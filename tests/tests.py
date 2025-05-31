import time
import os
from dataclasses import dataclass
from typing import List, Tuple
import json
from datetime import datetime
import matplotlib.pyplot as plt
import random

@dataclass
class TestResult:
    b_value: int
    s_value: int
    time: float
    algorithm: str  # 'mpqs', 'pmpqs', or 'dixon'

def test_configuration(numbers: List[int], b: int, s: int = 100000, algorithm: str = 'mpqs') -> float:
    """Test a specific configuration and return average execution time"""
    total_time = 0
    for n in numbers:
        start = time.time()
        os.system(f"../c/bin/factor -t {algorithm} -b {b} -s {s} -q -n {n}")
        end = time.time()
        total_time += end - start
    return total_time / len(numbers)

def binary_search_minimum(numbers: List[int], start: int, end: int, step: int, 
                         param: str = 'b', fixed_param: int = None, 
                         threshold: float = 0.001, algorithm: str = 'mpqs') -> Tuple[int, float]:
    """
    Binary search approach to find global minimum execution time
    param: 'b' for B-value search or 's' for S-value search
    fixed_param: fixed value for the other parameter (S when searching B, or B when searching S)
    """
    if end - start <= step:
        time_val = test_configuration(numbers, 
                                    b=start if param == 'b' else fixed_param,
                                    s=start if param == 's' else fixed_param,
                                    algorithm=algorithm)
        return start, time_val

    mid = start + (end - start) // 2
    points = [
        (mid - step, test_configuration(numbers,
                                      b=mid-step if param == 'b' else fixed_param,
                                      s=mid-step if param == 's' else fixed_param,
                                      algorithm=algorithm)),
        (mid, test_configuration(numbers,
                               b=mid if param == 'b' else fixed_param,
                               s=mid if param == 's' else fixed_param,
                               algorithm=algorithm)),
        (mid + step, test_configuration(numbers,
                                      b=mid+step if param == 'b' else fixed_param,
                                      s=mid+step if param == 's' else fixed_param,
                                      algorithm=algorithm))
    ]
    
    print(f"Testing {param}=({format(mid - step, ',')}, {format(mid, ',')}, {format(mid + step, ',')}) | times=({points[0][1]:.4f}s, {points[1][1]:.4f}s, {points[2][1]:.4f}s)")
    
    # Find the point with minimum time
    min_point = min(points, key=lambda x: x[1])
    
    # If difference between times is small enough, we're at a minimum
    if all(abs(p[1] - min_point[1]) < threshold for p in points):
        return min_point

    # Otherwise, continue searching in the direction of smaller times
    if points[0][1] < points[1][1]:
        return binary_search_minimum(numbers, start, mid - step, max(1, step // 2), param, fixed_param, threshold, algorithm)
    elif points[2][1] < points[1][1]:
        return binary_search_minimum(numbers, mid + step, end, max(1, step // 2), param, fixed_param, threshold, algorithm)
    else:
        # If we're in a valley, search with finer step
        if step > 100:
            return binary_search_minimum(numbers, mid - step, mid + step, max(1, step // 2), param, fixed_param, threshold, algorithm)
        return min_point

def get_interpolated_ranges(bit_count: int, algorithm: str) -> dict:
    """
    Interpolate good parameter ranges based on historical optimization results for specific algorithm
    """
    backup_file = os.path.join(os.path.dirname(__file__), "optimization_results.json")
    if not os.path.exists(backup_file):
        return None
        
    with open(backup_file, "r") as f:
        data = json.load(f)
    
    # Filter data for specific algorithm and sort by bit count
    algorithm_data = [entry for entry in data if entry.get("algorithm", "mpqs") == algorithm]
    sorted_data = sorted(algorithm_data, key=lambda x: x["bit_count"])
    
    if not sorted_data:
        return None
    
    # Find closest bit counts
    lower = None
    upper = None
    for entry in sorted_data:
        if entry["bit_count"] < bit_count:
            lower = entry
        if entry["bit_count"] > bit_count and upper is None:
            upper = entry

        # If we have exact match, use a range around that value
        if entry["bit_count"] == bit_count:
            b_val = entry["b_value"]
            return {
                'b_start': int(b_val * 0.8),
                'b_end': int(b_val * 1.2),
                'b_step': max(100, int(b_val * 0.05))
            }
    
    if lower is None or upper is None:
        return None

    # Linear interpolation
    bit_ratio = (bit_count - lower["bit_count"]) / (upper["bit_count"] - lower["bit_count"])
    interpolated_b = lower["b_value"] + bit_ratio * (upper["b_value"] - lower["b_value"])
    
    return {
        'b_start': int(interpolated_b * 0.8),
        'b_end': int(interpolated_b * 1.2),
        'b_step': max(100, int(interpolated_b * 0.05))
    }

def find_global_minimum(numbers: List[int], param_ranges: dict, algorithm: str = 'mpqs') -> TestResult:
    """Find global minimum by trying multiple starting points"""
    max_bits = len(bin(max(numbers))) - 2
    
    # Try to get interpolated ranges from historical data
    interpolated_ranges = get_interpolated_ranges(max_bits, algorithm)
    
    if interpolated_ranges is not None:
        param_ranges = interpolated_ranges
    elif param_ranges['b_start'] is None:
        # Fallback to default ranges based on number size and algorithm
        if algorithm == 'dixon':
            param_ranges['b_start'] = 12000
            param_ranges['b_end'] = 25000
            param_ranges['b_step'] = 2000
        else:  # mpqs or pmpqs
            param_ranges['b_start'] = max(1000, max_bits * 100)
            param_ranges['b_end'] = max(2000, max_bits * 200)
            param_ranges['b_step'] = max(100, max_bits * 10)

    # Try multiple starting points for B value first
    best_result = None
    num_starts = 3
    
    for i in range(num_starts):
        print("--")
        start_b = param_ranges['b_start'] + (i * (param_ranges['b_end'] - param_ranges['b_start'])) // (num_starts - 1)
        b_val, b_time = binary_search_minimum(numbers, start_b, param_ranges['b_end'], 
                                            param_ranges['b_step'], 'b', 1_000_000, algorithm=algorithm)
        
        # Now find optimal S value for this B (only for mpqs/pmpqs)
        if algorithm in ['qsieve', 'mpqs', 'pmpqs']:
            s_val, s_time = binary_search_minimum(numbers, 100_000, 10_000_000, 1_000_000, 's', b_val, algorithm=algorithm)
        else:  # dixon doesn't use sieving interval
            s_val, s_time = 0, b_time
        
        result = TestResult(b_val, s_val, s_time, algorithm)
        print(f"Search {i+1} found minimum at B={b_val}, S={s_val} with time={s_time:.4f}s\n--\n")
        
        if best_result is None or result.time < best_result.time:
            best_result = result
    
    return best_result

def main(numbers: List[int], factors: List[Tuple[int, int]], bitcount: int, algorithm: str = 'mpqs'):
    print(f"Optimizing parameters for {bitcount}-bit numbers using {algorithm}...")
    
    param_ranges = {
        'b_start': None,
        'b_end': None,
        'b_step': None
    }
    
    result = find_global_minimum(numbers, param_ranges, algorithm)
    
    # Save to backup file
    backup_data = {
        "timestamp": datetime.now().isoformat(),
        "bit_count": bitcount,
        "num_samples": len(numbers),
        "b_value": result.b_value,
        "avg_time": result.time,
        "algorithm": result.algorithm,
        "numbers": numbers,
        "factors": factors
    }

    # Only include s_value for non-Dixon algorithms
    if algorithm != 'dixon':
        backup_data["s_value"] = result.s_value
    
    backup_file = os.path.join(os.path.dirname(__file__), "optimization_results.json")
    
    # Read existing data or create empty list
    existing_data = []
    if os.path.exists(backup_file):
        with open(backup_file, "r") as f:
            try:
                existing_data = json.load(f)
            except json.JSONDecodeError:
                existing_data = []
    
    # Append new data and write back
    if not isinstance(existing_data, list):
        existing_data = [existing_data] if existing_data else []
    existing_data.append(backup_data)
    
    with open(backup_file, "w") as f:
        json.dump(existing_data, f, indent=4)
    print(f"Saved results to {backup_file}\n==\n")

    return result.b_value, result.s_value, result.time

def plot_optimization_results():
    """
    Read optimization results from JSON file, sort by bit size, and create plots
    showing average factoring time vs bit size for each algorithm, with uncertainty bars.
    """
    backup_file = os.path.join(os.path.dirname(__file__), "optimization_results.json")
    
    if not os.path.exists(backup_file):
        print("No optimization results file found")
        return
        
    with open(backup_file, "r") as f:
        data = json.load(f)
    
    # Group data by algorithm and bit count
    from collections import defaultdict
    grouped = defaultdict(lambda: defaultdict(list))
    for entry in data:
        algo = entry.get("algorithm", "mpqs")
        bit_count = entry["bit_count"]
        grouped[algo][bit_count].append(entry["avg_time"])
    
    plt.figure(figsize=(12, 8))
    
    for algo, bit_dict in grouped.items():
        bit_counts = []
        means = []
        lower_err = []
        upper_err = []
        for bit_count in sorted(bit_dict):
            times = bit_dict[bit_count]
            mean = sum(times) / len(times)
            min_time = min(times)
            max_time = max(times)
            bit_counts.append(bit_count)
            means.append(mean)
            lower_err.append(mean - min_time)
            upper_err.append(max_time - mean)
        # Error bars: asymmetric, from min to max
        plt.errorbar(
            bit_counts, means, 
            yerr=[lower_err, upper_err], 
            fmt='o-', capsize=5, label=algo.upper()
        )
    
    plt.yscale('log')
    plt.xlabel('Taille en bits')
    plt.ylabel('Temps moyen de factorisation (secondes)')
    plt.grid(True)
    plt.legend()
    
    plot_file = os.path.join(os.path.dirname(__file__), "optimization_plot.png")
    plt.savefig(plot_file)
    plt.close()
    
    print(f"Plot saved as {plot_file}")

def is_probably_prime(n: int, k: int = 5) -> bool:
    """Miller-Rabin primality test"""
    if n <= 3:
        return n > 1
    if n % 2 == 0:
        return False
    
    # Write n-1 as d * 2^r
    r = 0
    d = n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    
    # Witness loop
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True

def generate_prime(bits: int) -> int:
    """Generate a random prime number of specified bit length"""
    while True:
        # Generate random odd number of specified bit length
        n = random.getrandbits(bits) | (1 << (bits-1)) | 1
        if is_probably_prime(n, k=10):
            return n

def generate_semiprime(bit_length: int) -> Tuple[int, int, int]:
    """Generate a semiprime of specified bit length along with its prime factors"""
    # Split the bits roughly in half for the two primes
    p_bits = bit_length // 2
    q_bits = bit_length - p_bits
    
    while True:
        p = generate_prime(p_bits)
        q = generate_prime(q_bits)
        n = p * q
        actual_bits = len(bin(n)) - 2  # -2 for '0b' prefix
        
        # Ensure the product has exactly the desired bit length
        if actual_bits == bit_length:
            return n, p, q

def generate_test_cases(n: int, bit_length: int) -> List[int]:
    """Generate n semiprimes of specified bit length"""
    semiprimes = []
    factors = []
    for _ in range(n):
        semiprime, p, q = generate_semiprime(bit_length)
        semiprimes.append(semiprime)
        factors.append((p, q))
        #print(f"Generated {bit_length}-bit semiprime: {semiprime}")
        #print(f"Factors: {p} * {q}")
    return semiprimes, factors

if __name__ == "__main__":
    # Test different bit lengths with different algorithms
    
    #bits = [140]
    #algorithms = ['qsieve']
    #
    #numbers, factors = generate_test_cases(1, 250)
    #print(numbers, factors)
    #for bit_length in bits:
    #    numbers, factors = generate_test_cases(1, bit_length)  # Generate same numbers for fair comparison
    #    for algo in algorithms:
    #        main(numbers, factors, bit_length, algo)
    
    # Generate plot after all tests
    
    plot_optimization_results()
