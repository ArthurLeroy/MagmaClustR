"""
tactis2_worker.py — Entraîne et prédit TACTiS-2 pour UNE itération.

Appelé par run_benchmark.py avec :
    python tactis2_worker.py --iteration=1

Lit les données .rds (même format que le benchmark R MOMT/MT),
entraîne TACTiS-2 (phase 1 marginales DSF + phase 2 copule),
prédit sur chaque tâche test, sauvegarde modèle + prédictions.

Correspondance données R → TACTiS-2 :
  - datasets$train_data                  → data_train  (entraînement)
  - datasets$data_obs                    → data_obs    (identification Output_IDs)
  - datasets$pred_tasks_data[task_id]    → data_pred   (contexte par tâche test)
  - datasets$test_tasks_data[task_id]    → data_test   (vérité terrain par tâche)
  - datasets$test_task_ids               → liste des tâches à prédire
"""

import argparse
import os
import sys
import time
import random
import traceback

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import rdata
import matplotlib
import warnings # <-- AJOUTEZ CECI
matplotlib.use("Agg")

from tactis.model.tactis import TACTiS

# On rend rdata et les transformers muets sur les avertissements inoffensifs
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================
# Configuration
# ============================================================
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

SEED_BASE = 42

# Chemins serveur
BASE_XP1 = "/scratch/agrenoui/NeurIPS_experiments/Experience_1"
BASE_XP2 = "/scratch/agrenoui/NeurIPS_experiments/Experience_2"

N_OUT = 4
N_TRAIN = 30
N_PRED_MAX = 100

# DATASET_DIR = os.path.join(
#     BASE_XP1, f"n_out_{N_OUT}", "Datasets",
#     f"train_{N_TRAIN}_pred_{N_PRED_MAX}",
# )

DATASET_DIR = os.path.join(
    BASE_XP1, "Datasets", f"n_out_{N_OUT}",
    f"train_{N_TRAIN}_pred_{N_PRED_MAX}",
)
MODEL_DIR = os.path.join(BASE_XP2, "Models_TACTiS")
PRED_DIR = os.path.join(BASE_XP2, "Predictions_TACTiS")

# Hyperparamètres TACTiS-2
LR = 1e-3
WD = 1e-6
CLIP = 10.0
EPOCHS_P1 = 300
EPOCHS_P2 = 200
PATIENCE = 50
NUM_SAMPLES = 20
GAP_SIZE_MIN = 5
GAP_SIZE_MAX = 12


# ============================================================
# Lecture des .rds
# ============================================================
def load_rds_datasets(iteration):
    """Charge datasets_iter_{i}.rds et retourne la structure Python."""
    rds_path = os.path.join(DATASET_DIR, f"datasets_iter_{iteration}.rds")
    if not os.path.exists(rds_path):
        raise FileNotFoundError(f"Fichier introuvable : {rds_path}")
    return rdata.read_rds(rds_path)


def rdata_to_df(obj):
    """Convertit un objet rdata (data.frame R) en pandas DataFrame."""
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, dict):
        return pd.DataFrame(obj)
    return pd.DataFrame(dict(obj))


def ensure_types(df):
    """Force les types des colonnes standard."""
    for col in ["Input", "Output"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["Input_ID", "Output_ID", "Task_ID"]:
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df


def get_test_task_ids(parsed):
    """Extrait la liste des test_task_ids en strings."""
    ids = parsed["test_task_ids"]
    if hasattr(ids, "tolist"):
        ids = ids.tolist()
    return [str(x) for x in ids]


def get_named_list_df(named_list, key):
    """Récupère un DataFrame dans une liste nommée R (dict), cherche par str(key)."""
    if key in named_list:
        return ensure_types(rdata_to_df(named_list[key]))
    for k in named_list:
        if str(k) == str(key):
            return ensure_types(rdata_to_df(named_list[k]))
    return None


# ============================================================
# Fonctions utilitaires
# ============================================================
def normalize_per_variable(df, stats_map):
    """Normalise (x - mu) / sigma par Output_ID."""
    df_n = df.copy()
    for oid, s in stats_map.items():
        mu, sigma = s["mean"], s["std"]
        if sigma == 0 or pd.isna(sigma):
            sigma = 1.0
        sel = df_n["Output_ID"] == oid
        df_n.loc[sel, "Output"] = (df_n.loc[sel, "Output"] - mu) / sigma
    return df_n


def prepare_tensors(df, time_map, out_map):
    """Convertit un DataFrame multi-tâches en tenseurs (B, V, T)."""
    unique_tasks = df["Task_ID"].unique()
    task_map = {t: i for i, t in enumerate(unique_tasks)}
    B, V, T = len(unique_tasks), len(out_map), len(time_map)

    values = torch.zeros(B, V, T)
    mask = torch.zeros(B, V, T, dtype=torch.bool)

    b = df["Task_ID"].map(task_map).values
    v = df["Output_ID"].map(out_map).values
    t = df["Input"].map(time_map).values
    val = df["Output"].values

    valid = ~np.isnan(t) & ~np.isnan(v)
    b, v, t, val = b[valid], v[valid], t[valid].astype(int), val[valid]

    values[b, v, t] = torch.tensor(val, dtype=torch.float32)
    mask[b, v, t] = True
    return values, mask, unique_tasks


def prepare_tensors_single(df, time_map, out_map):
    """Tenseurs pour une seule tâche (B=1)."""
    V, T = len(out_map), len(time_map)
    values = torch.zeros(1, V, T)
    mask = torch.zeros(1, V, T, dtype=torch.bool)

    v = df["Output_ID"].map(out_map).values
    t = df["Input"].map(time_map).values
    val = df["Output"].values

    valid = ~np.isnan(t) & ~np.isnan(v)
    v_arr = v[valid].astype(int)
    t_arr = t[valid].astype(int)
    val_arr = val[valid]

    values[0, v_arr, t_arr] = torch.tensor(val_arr, dtype=torch.float32)
    mask[0, v_arr, t_arr] = True
    return values, mask


# ============================================================
# Encodage TACTiS
# ============================================================
def encode_flow(model, values, mask, time_idx):
    B, S, T = values.shape
    device = values.device
    se = model.flow_series_encoder(torch.arange(S, device=device))
    se = se[None, :, None, :].expand(B, -1, T, -1)
    mf = mask.float().unsqueeze(-1)
    vi = (values * mask.float()).unsqueeze(-1)
    x = torch.cat([vi, se, mf], dim=-1)
    x = model.flow_input_encoder(x)
    if model.input_encoding_normalization:
        x = x * model.flow_encoder_embedding_dim ** 0.5
    te = time_idx.unsqueeze(1).expand(-1, S, -1)
    x = model.flow_time_encoding(x, te.to(int))
    x = model.flow_encoder(x)
    return x


def encode_copula(model, values, mask, time_idx):
    B, S, T = values.shape
    device = values.device
    se = model.copula_series_encoder(torch.arange(S, device=device))
    se = se[None, :, None, :].expand(B, -1, T, -1)
    mf = mask.float().unsqueeze(-1)
    vi = (values * mask.float()).unsqueeze(-1)
    x = torch.cat([vi, se, mf], dim=-1)
    x = model.copula_input_encoder(x)
    if model.input_encoding_normalization:
        x = x * model.copula_encoder_embedding_dim ** 0.5
    te = time_idx.unsqueeze(1).expand(-1, S, -1)
    x = model.copula_time_encoding(x, te.to(int))
    x = model.copula_encoder(x)
    return x


def forward_loss(model, values, mask, time_idx):
    """Loss via decoder.loss(), normalisation 'both'."""
    flow_enc = encode_flow(model, values, mask, time_idx)
    cop_enc = (
        None if model.skip_copula
        else encode_copula(model, values, mask, time_idx)
    )
    _ = model.decoder.loss(
        flow_encoded=flow_enc, copula_encoded=cop_enc,
        mask=mask, true_value=values,
    )
    marg = model.decoder.marginal_logdet
    cop = model.decoder.copula_loss
    n_pred_total = (~mask[0]).sum().item()
    norm = max(1, n_pred_total)
    if isinstance(cop, torch.Tensor) and cop.dim() > 0:
        cop = cop / norm
    return marg / norm, cop


def make_train_mask(batch_size, n_vars, T, gap_start, gap_end, device):
    """True partout sauf le gap artificiel (= cibles)."""
    m = torch.ones(batch_size, n_vars, T, dtype=torch.bool, device=device)
    m[:, :, gap_start:gap_end] = False
    return m


def get_model_params(num_vars):
    """Paramètres d'architecture TACTiS-2."""
    return dict(
        num_series=num_vars,
        flow_series_embedding_dim=5,
        copula_series_embedding_dim=5,
        flow_input_encoder_layers=2,
        copula_input_encoder_layers=2,
        bagging_size=None,
        input_encoding_normalization=True,
        data_normalization="none",
        loss_normalization="both",
        positional_encoding={"dropout": 0.0},
        flow_encoder=dict(
            attention_layers=2, attention_heads=1,
            attention_dim=16, attention_feedforward_dim=16, dropout=0.0,
        ),
        copula_encoder=dict(
            attention_layers=2, attention_heads=1,
            attention_dim=16, attention_feedforward_dim=16, dropout=0.0,
        ),
        copula_decoder=dict(
            min_u=0.05, max_u=0.95,
            attentional_copula=dict(
                attention_heads=3, attention_layers=1,
                attention_dim=8, mlp_layers=2, mlp_dim=48,
                resolution=20, dropout=0.0,
            ),
            dsf_marginal=dict(
                mlp_layers=2, mlp_dim=48,
                flow_layers=2, flow_hid_dim=48,
            ),
        ),
        skip_copula=True,
        experiment_mode="interpolation",
    )


# ============================================================
# run_iteration : cœur du traitement
# ============================================================
def run_iteration(iteration):
    print(f"\n{'#' * 70}")
    print(f"#  ITERATION {iteration}")
    print(f"{'#' * 70}")

    # Reproductibilité
    seed_i = SEED_BASE + iteration
    random.seed(seed_i)
    torch.manual_seed(seed_i)
    np.random.seed(seed_i)

    # ========================================================
    # 1. Chargement du .rds
    # ========================================================
    parsed = load_rds_datasets(iteration)

    # --- data_train : données d'entraînement (toutes les tâches train) ---
    data_train = ensure_types(rdata_to_df(parsed["train_data"]))
    if "Cluster_ID" in data_train.columns:
        data_train = data_train.drop(columns=["Cluster_ID"])

    # --- data_obs : données d'observation (pour identifier Output_IDs) ---
    data_obs = ensure_types(rdata_to_df(parsed["data_obs"]))

    # --- pred_tasks_data : contexte observé par tâche test ---
    pred_tasks_data = parsed["pred_tasks_data"]

    # --- test_tasks_data : vérité terrain par tâche test ---
    test_tasks_data = parsed["test_tasks_data"]

    # --- test_task_ids : quelles tâches prédire ---
    test_task_ids = get_test_task_ids(parsed)

    print(f"  Tâches train : {data_train['Task_ID'].nunique()}")
    print(f"  Tâches test  : {len(test_task_ids)}")

    # ========================================================
    # 2. Identifier les variables (Output_ID)
    # ========================================================
    all_outputs = np.sort(np.unique(np.concatenate([
        data_obs["Output_ID"].values,
        data_train["Output_ID"].values,
    ])))
    output_map = {o: i for i, o in enumerate(all_outputs)}
    NUM_VARS = len(output_map)

    # ========================================================
    # 3. Construire la grille temporelle globale
    #    Union de : inputs train + tous les inputs pred/test
    # ========================================================
    all_inputs_list = [data_train["Input"].values]
    for tid_str in test_task_ids:
        df_pred = get_named_list_df(pred_tasks_data, tid_str)
        df_test = get_named_list_df(test_tasks_data, tid_str)
        if df_pred is not None:
            all_inputs_list.append(df_pred["Input"].values)
        if df_test is not None:
            all_inputs_list.append(df_test["Input"].values)

    local_inputs = np.sort(np.unique(np.concatenate(all_inputs_list)))
    local_time_map = {float(t): i for i, t in enumerate(local_inputs)}
    T_LOCAL = len(local_time_map)

    print(f"  Grille : T={T_LOCAL}, V={NUM_VARS}, Outputs={list(all_outputs)}")

    # ========================================================
    # 4. Normalisation (stats du train uniquement)
    # ========================================================
    train_stats = (
        data_train.groupby("Output_ID")["Output"]
        .agg(["mean", "std"]).to_dict("index")
    )
    data_train_norm = normalize_per_variable(data_train, train_stats)

    # ========================================================
    # 5. Tenseurs d'entraînement
    # ========================================================
    train_vals, train_mask, train_tasks = prepare_tensors(
        data_train_norm, local_time_map, output_map,
    )
    N_TRAIN = train_vals.shape[0]

    # ========================================================
    # 6. Gap range pour l'entraînement (plus grand bloc continu)
    # ========================================================
    all_vars_present = train_mask[0].all(dim=0)
    positions_both = torch.where(all_vars_present)[0]
    pb = positions_both.numpy()

    if len(pb) > 0:
        diffs = np.diff(pb)
        breaks = np.where(diffs > 1)[0]
        blocks = np.split(pb, breaks + 1)
        largest_block = max(blocks, key=len)
        GAP_RANGE_MIN = int(largest_block[0])
        GAP_RANGE_MAX = int(largest_block[-1])
    else:
        GAP_RANGE_MIN, GAP_RANGE_MAX = 0, T_LOCAL - 1

    VAL_GAP_SIZE = 8
    VAL_GAP_START = (GAP_RANGE_MIN + GAP_RANGE_MAX) // 2 - VAL_GAP_SIZE // 2
    VAL_GAP_END = VAL_GAP_START + VAL_GAP_SIZE

    time_idx = torch.arange(T_LOCAL, dtype=torch.int)

    # ========================================================
    # 7. Instanciation du modèle
    # ========================================================
    model = TACTiS(**get_model_params(NUM_VARS)).to(DEVICE)

    # ========================================================
    # PHASE 1 — Marginales (DSF)
    # ========================================================
    t_train_start = time.time()

    optimizer = optim.Adam(model.parameters(), lr=LR, weight_decay=WD)
    best_val_p1, best_ep_p1, best_st_p1 = float("inf"), -1, None

    for ep in range(EPOCHS_P1):
        model.train()
        optimizer.zero_grad()

        gs = random.randint(GAP_SIZE_MIN, GAP_SIZE_MAX)
        max_s = max(GAP_RANGE_MIN, GAP_RANGE_MAX - gs + 1)
        g_start = random.randint(GAP_RANGE_MIN, max_s)
        g_end = g_start + gs

        mask_tr = make_train_mask(N_TRAIN, NUM_VARS, T_LOCAL, g_start, g_end, DEVICE)
        ti = time_idx.unsqueeze(0).expand(N_TRAIN, -1).to(DEVICE)

        marg, _ = forward_loss(model, train_vals.to(DEVICE), mask_tr, ti)
        loss = (-marg).mean()
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=CLIP)
        optimizer.step()

        if (ep + 1) % 10 == 0:
            model.eval()
            with torch.no_grad():
                vm = make_train_mask(
                    N_TRAIN, NUM_VARS, T_LOCAL,
                    VAL_GAP_START, VAL_GAP_END, DEVICE,
                )
                vt = time_idx.unsqueeze(0).expand(N_TRAIN, -1).to(DEVICE)
                v_marg, _ = forward_loss(model, train_vals.to(DEVICE), vm, vt)
                v_loss = (-v_marg).mean().item()

            if v_loss < best_val_p1:
                best_val_p1, best_ep_p1 = v_loss, ep
                best_st_p1 = {
                    k: v.cpu().clone() for k, v in model.state_dict().items()
                }
            elif ep - best_ep_p1 > PATIENCE:
                print(f"  P1 early stop epoch {ep + 1} (best: {best_ep_p1 + 1})")
                break

    if best_st_p1:
        model.load_state_dict(best_st_p1)
    print(f"  P1 terminée. Best epoch {best_ep_p1 + 1}, val_loss={best_val_p1:.4f}")

    # ========================================================
    # PHASE 2 — Copule (AttentionalCopula)
    # ========================================================
    model.initialize_stage2()
    model.decoder.skip_copula = False
    model.to(DEVICE)

    freeze_keys = [
        "flow_series_encoder", "flow_time_encoding",
        "flow_input_encoder", "flow_encoder",
        "decoder.marginal",
    ]
    for name, p in model.named_parameters():
        if any(k in name for k in freeze_keys):
            p.requires_grad = False

    params_p2 = [p for p in model.parameters() if p.requires_grad]
    optimizer2 = optim.Adam(params_p2, lr=LR, weight_decay=WD)
    best_val_p2, best_ep_p2, best_st_p2 = float("inf"), -1, None

    for ep in range(EPOCHS_P2):
        model.train()
        optimizer2.zero_grad()

        gs = random.randint(GAP_SIZE_MIN, GAP_SIZE_MAX)
        max_s = max(GAP_RANGE_MIN, GAP_RANGE_MAX - gs + 1)
        g_start = random.randint(GAP_RANGE_MIN, max_s)
        g_end = g_start + gs

        mask_tr = make_train_mask(N_TRAIN, NUM_VARS, T_LOCAL, g_start, g_end, DEVICE)
        ti = time_idx.unsqueeze(0).expand(N_TRAIN, -1).to(DEVICE)

        _, cop = forward_loss(model, train_vals.to(DEVICE), mask_tr, ti)

        if isinstance(cop, torch.Tensor) and cop.dim() > 0:
            loss = cop.mean()
        else:
            loss = torch.tensor(0.0, device=DEVICE, requires_grad=True)

        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=CLIP)
        optimizer2.step()

        if (ep + 1) % 10 == 0:
            model.eval()
            with torch.no_grad():
                vm = make_train_mask(
                    N_TRAIN, NUM_VARS, T_LOCAL,
                    VAL_GAP_START, VAL_GAP_END, DEVICE,
                )
                vt = time_idx.unsqueeze(0).expand(N_TRAIN, -1).to(DEVICE)
                _, v_cop = forward_loss(model, train_vals.to(DEVICE), vm, vt)
                v_loss = (
                    v_cop.mean().item()
                    if isinstance(v_cop, torch.Tensor) and v_cop.dim() > 0
                    else 0.0
                )

            if v_loss < best_val_p2:
                best_val_p2, best_ep_p2 = v_loss, ep
                best_st_p2 = {
                    k: v.cpu().clone() for k, v in model.state_dict().items()
                }
            elif ep - best_ep_p2 > PATIENCE:
                print(f"  P2 early stop epoch {ep + 1} (best: {best_ep_p2 + 1})")
                break

    if best_st_p2:
        model.load_state_dict(best_st_p2)
    print(f"  P2 terminée. Best epoch {best_ep_p2 + 1}, val_loss={best_val_p2:.4f}")

    t_train_end = time.time()
    training_time = t_train_end - t_train_start
    print(f"  Temps d'entraînement : {training_time:.1f}s")

    # ========================================================
    # 8. Sauvegarde du modèle
    # ========================================================
    os.makedirs(MODEL_DIR, exist_ok=True)
    model_path = os.path.join(MODEL_DIR, f"model_iter_{iteration}.pth")
    torch.save({
        "model_state_dict": model.state_dict(),
        "model_params": get_model_params(NUM_VARS),
        "training_time_seconds": training_time,
        "best_epoch_p1": best_ep_p1 + 1,
        "best_val_loss_p1": best_val_p1,
        "best_epoch_p2": best_ep_p2 + 1,
        "best_val_loss_p2": best_val_p2,
        "train_stats": train_stats,
        "local_inputs": local_inputs,
        "output_map": output_map,
        "test_task_ids": test_task_ids,
    }, model_path)
    print(f"  Modèle sauvegardé : {model_path}")

    # ========================================================
    # 9. Inférence : boucle sur chaque tâche test
    #
    #    Pour chaque test_task_id :
    #      - pred_tasks_data[task_id] = données observées (contexte)
    #        → placées sur la grille avec mask=True
    #      - test_tasks_data[task_id] = vérité terrain
    #        → positions où extraire les prédictions
    #
    #    C'est l'équivalent de :
    #      R: pred_data_task → données observées pour la tâche
    #      R: test_grid_inputs → positions où prédire
    # ========================================================
    model.eval()
    t_pred_start = time.time()

    all_task_preds = []

    for tid_str in test_task_ids:
        # Récupérer data_pred (contexte) et data_test (vérité terrain)
        df_pred_raw = get_named_list_df(pred_tasks_data, tid_str)
        df_test_raw = get_named_list_df(test_tasks_data, tid_str)

        if df_pred_raw is None:
            print(f"  WARNING: pas de pred_tasks_data pour tâche {tid_str}, skip")
            continue
        if df_test_raw is None:
            print(f"  WARNING: pas de test_tasks_data pour tâche {tid_str}, skip")
            continue

        # Retirer Cluster_ID si présent
        if "Cluster_ID" in df_pred_raw.columns:
            df_pred_raw = df_pred_raw.drop(columns=["Cluster_ID"])

        # Normaliser les données de contexte avec les stats du train
        df_pred_norm = normalize_per_variable(df_pred_raw, train_stats)

        # Préparer les tenseurs : le contexte est masqué True (observé)
        inf_vals, inf_mask = prepare_tensors_single(
            df_pred_norm, local_time_map, output_map,
        )
        inf_vals = inf_vals.to(DEVICE)
        inf_mask = inf_mask.to(DEVICE)
        inf_ti = time_idx.unsqueeze(0).to(DEVICE)

        with torch.no_grad():
            flow_enc = encode_flow(model, inf_vals, inf_mask, inf_ti)
            cop_enc = encode_copula(model, inf_vals, inf_mask, inf_ti)
            samples = model.decoder.sample(
                num_samples=NUM_SAMPLES,
                flow_encoded=flow_enc,
                copula_encoded=cop_enc,
                mask=inf_mask,
                true_value=inf_vals,
            )

        # Extraire les prédictions aux positions de test pour chaque Output_ID
        for _, row in df_test_raw.iterrows():
            oid = str(row["Output_ID"])
            inp_val = float(row["Input"])
            true_val = float(row["Output"])

            if oid not in output_map:
                continue
            if float(inp_val) not in local_time_map:
                continue

            v_idx = output_map[oid]
            t_idx = local_time_map[float(inp_val)]

            # Stats de dénormalisation pour cet Output_ID
            mu = train_stats[oid]["mean"]
            sigma = train_stats[oid]["std"]
            if sigma == 0 or pd.isna(sigma):
                sigma = 1.0

            samp = samples[0, v_idx, t_idx, :].cpu()  # (NUM_SAMPLES,)
            samp_denorm = samp * sigma + mu

            all_task_preds.append({
                "Task_ID": tid_str,
                "Output_ID": oid,
                "Input": inp_val,
                "True": true_val,
                "Pred_mean": samp_denorm.mean().item(),
                "Pred_median": samp_denorm.median().item(),
                "Q05": samp_denorm.quantile(0.05).item(),
                "Q10": samp_denorm.quantile(0.10).item(),
                "Q25": samp_denorm.quantile(0.25).item(),
                "Q75": samp_denorm.quantile(0.75).item(),
                "Q90": samp_denorm.quantile(0.90).item(),
                "Q95": samp_denorm.quantile(0.95).item(),
            })

    t_pred_end = time.time()
    prediction_time = t_pred_end - t_pred_start
    print(f"  Temps de prédiction : {prediction_time:.2f}s")

    # ========================================================
    # 10. Sauvegarde des prédictions
    # ========================================================
    os.makedirs(PRED_DIR, exist_ok=True)
    preds_df = pd.DataFrame(all_task_preds)
    preds_df["training_time_s"] = training_time
    preds_df["prediction_time_s"] = prediction_time

    pred_path = os.path.join(PRED_DIR, f"predictions_iter_{iteration}.csv")
    preds_df.to_csv(pred_path, index=False)
    print(f"  Prédictions sauvegardées : {pred_path}")

    # Métriques résumées
    if len(preds_df) > 0:
        valid = preds_df.dropna(subset=["True", "Pred_mean"])
        if len(valid) > 0:
            rmse = float(np.sqrt(((valid["Pred_mean"] - valid["True"]) ** 2).mean()))
            mae = float(np.abs(valid["Pred_mean"] - valid["True"]).mean())
            print(f"  RMSE={rmse:.4f}, MAE={mae:.4f}")

    # Nettoyage
    del model, optimizer, optimizer2, samples, flow_enc, cop_enc
    del train_vals, train_mask
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    print(f"  Iteration {iteration} terminée avec succès.")
    return 0


# ============================================================
# Point d'entrée
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TACTiS-2 worker — single iteration")
    parser.add_argument("--iteration", type=int, required=True,
                        help="Numéro de l'itération (1..100)")
    args = parser.parse_args()

    try:
        sys.exit(run_iteration(args.iteration))
    except Exception as e:
        print(f"\n*** ERREUR iteration {args.iteration} ***")
        print(f"{type(e).__name__}: {e}")
        traceback.print_exc()
        sys.exit(1)
