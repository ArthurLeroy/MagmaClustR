"""
run_benchmark.py — Orchestrateur parallèle pour le benchmark TACTiS-2.

Lance les 100 itérations sur 16 workers en parallèle avec du
work-stealing : dès qu'un worker finit, il prend l'itération suivante.

Usage :
    python run_benchmark.py

Chaque itération est lancée comme un sous-processus indépendant
(tactis2_worker.py) pour garantir l'isolation mémoire.
"""

import os
import sys
import time
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

# ============================================================
# Configuration
# ============================================================
N_ITERATIONS = 100
N_WORKERS = 32

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WORKER_SCRIPT = os.path.join(SCRIPT_DIR, "tactis2_worker.py")

BASE_XP2 = "/scratch/agrenoui/NeurIPS_experiments/Experience_2"
LOG_DIR = os.path.join(BASE_XP2, "logs_TACTiS")
os.makedirs(LOG_DIR, exist_ok=True)


def run_single_iteration(iteration):
    """Lance tactis2_worker.py pour une itération et capture le log."""
    log_file = os.path.join(LOG_DIR, f"iter_{iteration}.log")
    t_start = time.time()

    try:
        with open(log_file, "w") as f_log:
            result = subprocess.run(
                [sys.executable, "-u", WORKER_SCRIPT, f"--iteration={iteration}"],
                stdout=f_log,
                stderr=subprocess.STDOUT,
                timeout=18000,  # 5h max par itération
            )
        elapsed = time.time() - t_start
        status = "OK" if result.returncode == 0 else f"ERREUR (code {result.returncode})"
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t_start
        status = "TIMEOUT"
    except Exception as e:
        elapsed = time.time() - t_start
        status = f"EXCEPTION: {e}"

    return iteration, status, elapsed


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print(f"{'=' * 70}")
    print(f"  Benchmark TACTiS-2 — {N_ITERATIONS} itérations, {N_WORKERS} workers")
    print(f"  Logs : {LOG_DIR}")
    print(f"  Début : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 70}\n")

    t_global_start = time.time()
    results = []

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        # Soumettre toutes les itérations d'un coup.
        # L'executor gère le work-stealing : dès qu'un worker est libre,
        # il prend le prochain job dans la queue.
        futures = {
            executor.submit(run_single_iteration, i): i
            for i in range(1, N_ITERATIONS + 1)
        }

        for future in as_completed(futures):
            iteration, status, elapsed = future.result()
            results.append((iteration, status, elapsed))
            n_done = len(results)
            print(
                f"  [{n_done:3d}/{N_ITERATIONS}] Iter {iteration:3d} : "
                f"{status:<20s} ({elapsed:.1f}s)"
            )

    # ============================================================
    # Récapitulatif
    # ============================================================
    t_global = time.time() - t_global_start
    results.sort(key=lambda x: x[0])

    n_ok = sum(1 for _, s, _ in results if s == "OK")
    n_err = sum(1 for _, s, _ in results if s != "OK")

    # Sauvegarder le résumé
    import csv
    summary_path = os.path.join(BASE_XP2, "benchmark_summary_TACTiS2.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["iteration", "status", "elapsed_s"])
        for it, st, el in results:
            writer.writerow([it, st, f"{el:.1f}"])

    print(f"\n{'=' * 70}")
    print(f"  Benchmark terminé.")
    print(f"  Réussies : {n_ok}/{N_ITERATIONS}")
    print(f"  Erreurs  : {n_err}/{N_ITERATIONS}")
    print(f"  Durée totale : {t_global / 3600:.2f}h")
    print(f"  Récapitulatif : {summary_path}")
    print(f"{'=' * 70}")
