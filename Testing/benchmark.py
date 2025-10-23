import subprocess # for running subprocesses and capturing output
import sys
import time
import statistics as stats
from pathlib import Path
import argparse
import csv
import notebook_utils

PY = sys.executable
# use the same interpreter to run scripts, not just a random python

ROOT = Path(__file__).resolve().parent  # important for subprocess, 

def human(s: float) -> str: # in order to make output human-readable with µs and ms
    if s < 1e-3: return f"{s*1e6:.1f} µs"
    if s < 1:   return f"{s*1e3:.1f} ms"
    return f"{s:.3f} s"

def run_once(script_path: Path, *args, timeout=500) -> float: #needed for running the script once and measuring time, because of the subprocess
    t0 = time.perf_counter() # for measuring elapsed time, it makes precise measurements
    process = subprocess.run([PY, str(script_path), *map(str, args)],
                  cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=timeout)
    dt = time.perf_counter() - t0
    # actually the timing is done here now capturing the errors
    if process.returncode != 0:
        err = process.stderr.decode(errors="replace")
        raise RuntimeError(f"{script_path.name} error:\n{err}")
    return dt

def benchmark_script(script_rel: str, rep: int = 5, warmup: int = 1, timeout: int = 600, *args): # timeout 12 minutes

    script = ROOT / script_rel

    # if the script does not exist, skip it
    if not script.exists():
        print("script does not exist\n")
        return None

    # warmup runs
    for _ in range(warmup):
        try:
            run_once(script, *args, timeout=timeout)
        except Exception as e:
            print("[warmup error]")
            break


    times = []
    t_total0 = time.perf_counter() #starts
    for i in range(rep):
        dt = run_once(script, *args, timeout=timeout)
        times.append(dt)
        #print(f"Run {i+1}/{repeats} finished in {human(dt)}")
    total = time.perf_counter() - t_total0

    if not times:
        print(f"[warn] no measurement: {script_rel}")
        return None

    best = min(times)
    avg  = stats.mean(times)
    std  = stats.pstdev(times) if len(times) > 1 else 0.0
    rps  = len(times) / total
    per1000 = 1000.0 / rps

    print(f"\nResults for {script_rel}:")
    print(f"  Best time: {human(best)}")
    print(f"  Average time: {human(avg)}")
    print(f"  Std deviation: {human(std)}")
    print(f"  Total time: {human(total)}")
    print(f"  Runs per second: {rps:.2f}")
    print(f"  Time for 1000 runs: {human(per1000)}\n")


    return {
        "script": script_rel,
        "runs": len(times),
        "warmup": warmup,
        "best_s": best,
        "avg_s": avg,
        "stdev_s": std,
        "total_s": total,
        "runs_per_sec": rps,
        "per_1000_s": per1000,
    }

def main():
    p = argparse.ArgumentParser() # for parsing command-line arguments
    p.add_argument("-rep", type=int, default=100, help="how many times")
    p.add_argument("-warmup", type=int, default=1, help="warmup runs")
    p.add_argument("-timeout", type=int, default=600, help="time out in seconds")
    p.add_argument("-csv", type=str, default="bench_results.csv", help="CSV output file")
    p.add_argument("scripts", nargs="*", default=["potentials/main.py", "potentials/main_gravitation.py", "potentials/mie_potential.py"],
                   help="measured scripts")
    args = p.parse_args()

    rows = []
    for s in args.scripts:
        res = benchmark_script(s, rep=args.rep, warmup=args.warmup, timeout=args.timeout)
        if res: rows.append(res)

    if rows:
        out = Path(args.csv)
        with out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\nSaved CSV: {out.resolve()}")
        notebook_utils.make_notebook(out)
    else:
        print("\n[info] yazılacak sonuç yok.")

if __name__ == "__main__":
    main()
