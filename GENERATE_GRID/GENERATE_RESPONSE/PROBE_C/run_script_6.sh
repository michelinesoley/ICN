#!/bin/bash
#SBATCH -J run_script_6
#SBATCH --partition=pi_esi
#SBATCH --time 120:00:00
#SBATCH --ntasks=12 --nodes=1
#SBATCH --mem-per-cpu=5000
#SBATCH -o run_script_6.err

### Load MPI module
#module load MPI/OpenMPI/2.1.1-intel15

### Load QCHEM module
module load Q-Chem/5.2-openmpi

export QCSCRATCH=/gpfs/loomis/scratch60/fas/batista/pev4/qchem-scratch/
 
### Run program

echo 'Starting program on'
date
echo ''

qchem -nt $SLURM_NTASKS qchem_2205.in qchem_2205.out
qchem -nt $SLURM_NTASKS qchem_2206.in qchem_2206.out
qchem -nt $SLURM_NTASKS qchem_2207.in qchem_2207.out
qchem -nt $SLURM_NTASKS qchem_2208.in qchem_2208.out
qchem -nt $SLURM_NTASKS qchem_2209.in qchem_2209.out
qchem -nt $SLURM_NTASKS qchem_2210.in qchem_2210.out
qchem -nt $SLURM_NTASKS qchem_2211.in qchem_2211.out
qchem -nt $SLURM_NTASKS qchem_2212.in qchem_2212.out
qchem -nt $SLURM_NTASKS qchem_2213.in qchem_2213.out
qchem -nt $SLURM_NTASKS qchem_2214.in qchem_2214.out
qchem -nt $SLURM_NTASKS qchem_2215.in qchem_2215.out
qchem -nt $SLURM_NTASKS qchem_2216.in qchem_2216.out
qchem -nt $SLURM_NTASKS qchem_2217.in qchem_2217.out
qchem -nt $SLURM_NTASKS qchem_2218.in qchem_2218.out
qchem -nt $SLURM_NTASKS qchem_2219.in qchem_2219.out
qchem -nt $SLURM_NTASKS qchem_2220.in qchem_2220.out
qchem -nt $SLURM_NTASKS qchem_2221.in qchem_2221.out
qchem -nt $SLURM_NTASKS qchem_2222.in qchem_2222.out
qchem -nt $SLURM_NTASKS qchem_2223.in qchem_2223.out
qchem -nt $SLURM_NTASKS qchem_2224.in qchem_2224.out
qchem -nt $SLURM_NTASKS qchem_2225.in qchem_2225.out
qchem -nt $SLURM_NTASKS qchem_2226.in qchem_2226.out
qchem -nt $SLURM_NTASKS qchem_2227.in qchem_2227.out
qchem -nt $SLURM_NTASKS qchem_2228.in qchem_2228.out
qchem -nt $SLURM_NTASKS qchem_2229.in qchem_2229.out
qchem -nt $SLURM_NTASKS qchem_2230.in qchem_2230.out
qchem -nt $SLURM_NTASKS qchem_2231.in qchem_2231.out
qchem -nt $SLURM_NTASKS qchem_2232.in qchem_2232.out
qchem -nt $SLURM_NTASKS qchem_2233.in qchem_2233.out
qchem -nt $SLURM_NTASKS qchem_2234.in qchem_2234.out
qchem -nt $SLURM_NTASKS qchem_2235.in qchem_2235.out
qchem -nt $SLURM_NTASKS qchem_2236.in qchem_2236.out
qchem -nt $SLURM_NTASKS qchem_2237.in qchem_2237.out
qchem -nt $SLURM_NTASKS qchem_2238.in qchem_2238.out
qchem -nt $SLURM_NTASKS qchem_2239.in qchem_2239.out
qchem -nt $SLURM_NTASKS qchem_2240.in qchem_2240.out
qchem -nt $SLURM_NTASKS qchem_2241.in qchem_2241.out
qchem -nt $SLURM_NTASKS qchem_2242.in qchem_2242.out
qchem -nt $SLURM_NTASKS qchem_2243.in qchem_2243.out
qchem -nt $SLURM_NTASKS qchem_2244.in qchem_2244.out
qchem -nt $SLURM_NTASKS qchem_2245.in qchem_2245.out
qchem -nt $SLURM_NTASKS qchem_2246.in qchem_2246.out
qchem -nt $SLURM_NTASKS qchem_2247.in qchem_2247.out
qchem -nt $SLURM_NTASKS qchem_2248.in qchem_2248.out
qchem -nt $SLURM_NTASKS qchem_2249.in qchem_2249.out
qchem -nt $SLURM_NTASKS qchem_2250.in qchem_2250.out
qchem -nt $SLURM_NTASKS qchem_2251.in qchem_2251.out
qchem -nt $SLURM_NTASKS qchem_2252.in qchem_2252.out
qchem -nt $SLURM_NTASKS qchem_2253.in qchem_2253.out
qchem -nt $SLURM_NTASKS qchem_2254.in qchem_2254.out
qchem -nt $SLURM_NTASKS qchem_2255.in qchem_2255.out
qchem -nt $SLURM_NTASKS qchem_2256.in qchem_2256.out
qchem -nt $SLURM_NTASKS qchem_2257.in qchem_2257.out
qchem -nt $SLURM_NTASKS qchem_2258.in qchem_2258.out
qchem -nt $SLURM_NTASKS qchem_2259.in qchem_2259.out
qchem -nt $SLURM_NTASKS qchem_2260.in qchem_2260.out
qchem -nt $SLURM_NTASKS qchem_2261.in qchem_2261.out
qchem -nt $SLURM_NTASKS qchem_2262.in qchem_2262.out
qchem -nt $SLURM_NTASKS qchem_2263.in qchem_2263.out
qchem -nt $SLURM_NTASKS qchem_2264.in qchem_2264.out
qchem -nt $SLURM_NTASKS qchem_2265.in qchem_2265.out
qchem -nt $SLURM_NTASKS qchem_2266.in qchem_2266.out
qchem -nt $SLURM_NTASKS qchem_2267.in qchem_2267.out
qchem -nt $SLURM_NTASKS qchem_2268.in qchem_2268.out
qchem -nt $SLURM_NTASKS qchem_2269.in qchem_2269.out
qchem -nt $SLURM_NTASKS qchem_2270.in qchem_2270.out
qchem -nt $SLURM_NTASKS qchem_2271.in qchem_2271.out
qchem -nt $SLURM_NTASKS qchem_2272.in qchem_2272.out
qchem -nt $SLURM_NTASKS qchem_2273.in qchem_2273.out
qchem -nt $SLURM_NTASKS qchem_2274.in qchem_2274.out
qchem -nt $SLURM_NTASKS qchem_2275.in qchem_2275.out
qchem -nt $SLURM_NTASKS qchem_2276.in qchem_2276.out
qchem -nt $SLURM_NTASKS qchem_2277.in qchem_2277.out
qchem -nt $SLURM_NTASKS qchem_2278.in qchem_2278.out
qchem -nt $SLURM_NTASKS qchem_2279.in qchem_2279.out
qchem -nt $SLURM_NTASKS qchem_2280.in qchem_2280.out
qchem -nt $SLURM_NTASKS qchem_2281.in qchem_2281.out
qchem -nt $SLURM_NTASKS qchem_2282.in qchem_2282.out
qchem -nt $SLURM_NTASKS qchem_2283.in qchem_2283.out
qchem -nt $SLURM_NTASKS qchem_2284.in qchem_2284.out
qchem -nt $SLURM_NTASKS qchem_2285.in qchem_2285.out
qchem -nt $SLURM_NTASKS qchem_2286.in qchem_2286.out
qchem -nt $SLURM_NTASKS qchem_2287.in qchem_2287.out
qchem -nt $SLURM_NTASKS qchem_2288.in qchem_2288.out
qchem -nt $SLURM_NTASKS qchem_2289.in qchem_2289.out
qchem -nt $SLURM_NTASKS qchem_2290.in qchem_2290.out
qchem -nt $SLURM_NTASKS qchem_2291.in qchem_2291.out
qchem -nt $SLURM_NTASKS qchem_2292.in qchem_2292.out
qchem -nt $SLURM_NTASKS qchem_2293.in qchem_2293.out
qchem -nt $SLURM_NTASKS qchem_2294.in qchem_2294.out
qchem -nt $SLURM_NTASKS qchem_2295.in qchem_2295.out
qchem -nt $SLURM_NTASKS qchem_2296.in qchem_2296.out
qchem -nt $SLURM_NTASKS qchem_2297.in qchem_2297.out
qchem -nt $SLURM_NTASKS qchem_2298.in qchem_2298.out
qchem -nt $SLURM_NTASKS qchem_2299.in qchem_2299.out
qchem -nt $SLURM_NTASKS qchem_2300.in qchem_2300.out
qchem -nt $SLURM_NTASKS qchem_2301.in qchem_2301.out
qchem -nt $SLURM_NTASKS qchem_2302.in qchem_2302.out
qchem -nt $SLURM_NTASKS qchem_2303.in qchem_2303.out
qchem -nt $SLURM_NTASKS qchem_2304.in qchem_2304.out
qchem -nt $SLURM_NTASKS qchem_2305.in qchem_2305.out
qchem -nt $SLURM_NTASKS qchem_2306.in qchem_2306.out
qchem -nt $SLURM_NTASKS qchem_2307.in qchem_2307.out
qchem -nt $SLURM_NTASKS qchem_2308.in qchem_2308.out
qchem -nt $SLURM_NTASKS qchem_2309.in qchem_2309.out
qchem -nt $SLURM_NTASKS qchem_2310.in qchem_2310.out
qchem -nt $SLURM_NTASKS qchem_2311.in qchem_2311.out
qchem -nt $SLURM_NTASKS qchem_2312.in qchem_2312.out
qchem -nt $SLURM_NTASKS qchem_2313.in qchem_2313.out
qchem -nt $SLURM_NTASKS qchem_2314.in qchem_2314.out
qchem -nt $SLURM_NTASKS qchem_2315.in qchem_2315.out
qchem -nt $SLURM_NTASKS qchem_2316.in qchem_2316.out
qchem -nt $SLURM_NTASKS qchem_2317.in qchem_2317.out
qchem -nt $SLURM_NTASKS qchem_2318.in qchem_2318.out
qchem -nt $SLURM_NTASKS qchem_2319.in qchem_2319.out
qchem -nt $SLURM_NTASKS qchem_2320.in qchem_2320.out
qchem -nt $SLURM_NTASKS qchem_2321.in qchem_2321.out
qchem -nt $SLURM_NTASKS qchem_2322.in qchem_2322.out
qchem -nt $SLURM_NTASKS qchem_2323.in qchem_2323.out
qchem -nt $SLURM_NTASKS qchem_2324.in qchem_2324.out
qchem -nt $SLURM_NTASKS qchem_2325.in qchem_2325.out
qchem -nt $SLURM_NTASKS qchem_2326.in qchem_2326.out
qchem -nt $SLURM_NTASKS qchem_2327.in qchem_2327.out
qchem -nt $SLURM_NTASKS qchem_2328.in qchem_2328.out
qchem -nt $SLURM_NTASKS qchem_2329.in qchem_2329.out
qchem -nt $SLURM_NTASKS qchem_2330.in qchem_2330.out
qchem -nt $SLURM_NTASKS qchem_2331.in qchem_2331.out
qchem -nt $SLURM_NTASKS qchem_2332.in qchem_2332.out
qchem -nt $SLURM_NTASKS qchem_2333.in qchem_2333.out
qchem -nt $SLURM_NTASKS qchem_2334.in qchem_2334.out
qchem -nt $SLURM_NTASKS qchem_2335.in qchem_2335.out
qchem -nt $SLURM_NTASKS qchem_2336.in qchem_2336.out
qchem -nt $SLURM_NTASKS qchem_2337.in qchem_2337.out
qchem -nt $SLURM_NTASKS qchem_2338.in qchem_2338.out
qchem -nt $SLURM_NTASKS qchem_2339.in qchem_2339.out
qchem -nt $SLURM_NTASKS qchem_2340.in qchem_2340.out
qchem -nt $SLURM_NTASKS qchem_2341.in qchem_2341.out
qchem -nt $SLURM_NTASKS qchem_2342.in qchem_2342.out
qchem -nt $SLURM_NTASKS qchem_2343.in qchem_2343.out
qchem -nt $SLURM_NTASKS qchem_2344.in qchem_2344.out
qchem -nt $SLURM_NTASKS qchem_2345.in qchem_2345.out
qchem -nt $SLURM_NTASKS qchem_2346.in qchem_2346.out
qchem -nt $SLURM_NTASKS qchem_2347.in qchem_2347.out
qchem -nt $SLURM_NTASKS qchem_2348.in qchem_2348.out
qchem -nt $SLURM_NTASKS qchem_2349.in qchem_2349.out
qchem -nt $SLURM_NTASKS qchem_2350.in qchem_2350.out
qchem -nt $SLURM_NTASKS qchem_2351.in qchem_2351.out
qchem -nt $SLURM_NTASKS qchem_2352.in qchem_2352.out
qchem -nt $SLURM_NTASKS qchem_2353.in qchem_2353.out
qchem -nt $SLURM_NTASKS qchem_2354.in qchem_2354.out
qchem -nt $SLURM_NTASKS qchem_2355.in qchem_2355.out
qchem -nt $SLURM_NTASKS qchem_2356.in qchem_2356.out
qchem -nt $SLURM_NTASKS qchem_2357.in qchem_2357.out
qchem -nt $SLURM_NTASKS qchem_2358.in qchem_2358.out
qchem -nt $SLURM_NTASKS qchem_2359.in qchem_2359.out
qchem -nt $SLURM_NTASKS qchem_2360.in qchem_2360.out
qchem -nt $SLURM_NTASKS qchem_2361.in qchem_2361.out
qchem -nt $SLURM_NTASKS qchem_2362.in qchem_2362.out
qchem -nt $SLURM_NTASKS qchem_2363.in qchem_2363.out
qchem -nt $SLURM_NTASKS qchem_2364.in qchem_2364.out
qchem -nt $SLURM_NTASKS qchem_2365.in qchem_2365.out
qchem -nt $SLURM_NTASKS qchem_2366.in qchem_2366.out
qchem -nt $SLURM_NTASKS qchem_2367.in qchem_2367.out
qchem -nt $SLURM_NTASKS qchem_2368.in qchem_2368.out
qchem -nt $SLURM_NTASKS qchem_2369.in qchem_2369.out
qchem -nt $SLURM_NTASKS qchem_2370.in qchem_2370.out
qchem -nt $SLURM_NTASKS qchem_2371.in qchem_2371.out
qchem -nt $SLURM_NTASKS qchem_2372.in qchem_2372.out
qchem -nt $SLURM_NTASKS qchem_2373.in qchem_2373.out
qchem -nt $SLURM_NTASKS qchem_2374.in qchem_2374.out
qchem -nt $SLURM_NTASKS qchem_2375.in qchem_2375.out
qchem -nt $SLURM_NTASKS qchem_2376.in qchem_2376.out
qchem -nt $SLURM_NTASKS qchem_2377.in qchem_2377.out
qchem -nt $SLURM_NTASKS qchem_2378.in qchem_2378.out
qchem -nt $SLURM_NTASKS qchem_2379.in qchem_2379.out
qchem -nt $SLURM_NTASKS qchem_2380.in qchem_2380.out
qchem -nt $SLURM_NTASKS qchem_2381.in qchem_2381.out
qchem -nt $SLURM_NTASKS qchem_2382.in qchem_2382.out
qchem -nt $SLURM_NTASKS qchem_2383.in qchem_2383.out
qchem -nt $SLURM_NTASKS qchem_2384.in qchem_2384.out
qchem -nt $SLURM_NTASKS qchem_2385.in qchem_2385.out
qchem -nt $SLURM_NTASKS qchem_2386.in qchem_2386.out
qchem -nt $SLURM_NTASKS qchem_2387.in qchem_2387.out
qchem -nt $SLURM_NTASKS qchem_2388.in qchem_2388.out
qchem -nt $SLURM_NTASKS qchem_2389.in qchem_2389.out
qchem -nt $SLURM_NTASKS qchem_2390.in qchem_2390.out
qchem -nt $SLURM_NTASKS qchem_2391.in qchem_2391.out
qchem -nt $SLURM_NTASKS qchem_2392.in qchem_2392.out
qchem -nt $SLURM_NTASKS qchem_2393.in qchem_2393.out
qchem -nt $SLURM_NTASKS qchem_2394.in qchem_2394.out
qchem -nt $SLURM_NTASKS qchem_2395.in qchem_2395.out
qchem -nt $SLURM_NTASKS qchem_2396.in qchem_2396.out
qchem -nt $SLURM_NTASKS qchem_2397.in qchem_2397.out
qchem -nt $SLURM_NTASKS qchem_2398.in qchem_2398.out
qchem -nt $SLURM_NTASKS qchem_2399.in qchem_2399.out
qchem -nt $SLURM_NTASKS qchem_2400.in qchem_2400.out
qchem -nt $SLURM_NTASKS qchem_2401.in qchem_2401.out
qchem -nt $SLURM_NTASKS qchem_2402.in qchem_2402.out
qchem -nt $SLURM_NTASKS qchem_2403.in qchem_2403.out
qchem -nt $SLURM_NTASKS qchem_2404.in qchem_2404.out
qchem -nt $SLURM_NTASKS qchem_2405.in qchem_2405.out
qchem -nt $SLURM_NTASKS qchem_2406.in qchem_2406.out
qchem -nt $SLURM_NTASKS qchem_2407.in qchem_2407.out
qchem -nt $SLURM_NTASKS qchem_2408.in qchem_2408.out
qchem -nt $SLURM_NTASKS qchem_2409.in qchem_2409.out
qchem -nt $SLURM_NTASKS qchem_2410.in qchem_2410.out
qchem -nt $SLURM_NTASKS qchem_2411.in qchem_2411.out
qchem -nt $SLURM_NTASKS qchem_2412.in qchem_2412.out
qchem -nt $SLURM_NTASKS qchem_2413.in qchem_2413.out
qchem -nt $SLURM_NTASKS qchem_2414.in qchem_2414.out
qchem -nt $SLURM_NTASKS qchem_2415.in qchem_2415.out
qchem -nt $SLURM_NTASKS qchem_2416.in qchem_2416.out
qchem -nt $SLURM_NTASKS qchem_2417.in qchem_2417.out
qchem -nt $SLURM_NTASKS qchem_2418.in qchem_2418.out
qchem -nt $SLURM_NTASKS qchem_2419.in qchem_2419.out
qchem -nt $SLURM_NTASKS qchem_2420.in qchem_2420.out
qchem -nt $SLURM_NTASKS qchem_2421.in qchem_2421.out
qchem -nt $SLURM_NTASKS qchem_2422.in qchem_2422.out
qchem -nt $SLURM_NTASKS qchem_2423.in qchem_2423.out
qchem -nt $SLURM_NTASKS qchem_2424.in qchem_2424.out
qchem -nt $SLURM_NTASKS qchem_2425.in qchem_2425.out
qchem -nt $SLURM_NTASKS qchem_2426.in qchem_2426.out
qchem -nt $SLURM_NTASKS qchem_2427.in qchem_2427.out
qchem -nt $SLURM_NTASKS qchem_2428.in qchem_2428.out
qchem -nt $SLURM_NTASKS qchem_2429.in qchem_2429.out
qchem -nt $SLURM_NTASKS qchem_2430.in qchem_2430.out
qchem -nt $SLURM_NTASKS qchem_2431.in qchem_2431.out
qchem -nt $SLURM_NTASKS qchem_2432.in qchem_2432.out
qchem -nt $SLURM_NTASKS qchem_2433.in qchem_2433.out
qchem -nt $SLURM_NTASKS qchem_2434.in qchem_2434.out
qchem -nt $SLURM_NTASKS qchem_2435.in qchem_2435.out
qchem -nt $SLURM_NTASKS qchem_2436.in qchem_2436.out
qchem -nt $SLURM_NTASKS qchem_2437.in qchem_2437.out
qchem -nt $SLURM_NTASKS qchem_2438.in qchem_2438.out
qchem -nt $SLURM_NTASKS qchem_2439.in qchem_2439.out
qchem -nt $SLURM_NTASKS qchem_2440.in qchem_2440.out
qchem -nt $SLURM_NTASKS qchem_2441.in qchem_2441.out
qchem -nt $SLURM_NTASKS qchem_2442.in qchem_2442.out
qchem -nt $SLURM_NTASKS qchem_2443.in qchem_2443.out
qchem -nt $SLURM_NTASKS qchem_2444.in qchem_2444.out
qchem -nt $SLURM_NTASKS qchem_2445.in qchem_2445.out
qchem -nt $SLURM_NTASKS qchem_2446.in qchem_2446.out
qchem -nt $SLURM_NTASKS qchem_2447.in qchem_2447.out
qchem -nt $SLURM_NTASKS qchem_2448.in qchem_2448.out
qchem -nt $SLURM_NTASKS qchem_2449.in qchem_2449.out
qchem -nt $SLURM_NTASKS qchem_2450.in qchem_2450.out
qchem -nt $SLURM_NTASKS qchem_2451.in qchem_2451.out
qchem -nt $SLURM_NTASKS qchem_2452.in qchem_2452.out
qchem -nt $SLURM_NTASKS qchem_2453.in qchem_2453.out
qchem -nt $SLURM_NTASKS qchem_2454.in qchem_2454.out
qchem -nt $SLURM_NTASKS qchem_2455.in qchem_2455.out
qchem -nt $SLURM_NTASKS qchem_2456.in qchem_2456.out
qchem -nt $SLURM_NTASKS qchem_2457.in qchem_2457.out
qchem -nt $SLURM_NTASKS qchem_2458.in qchem_2458.out
qchem -nt $SLURM_NTASKS qchem_2459.in qchem_2459.out
qchem -nt $SLURM_NTASKS qchem_2460.in qchem_2460.out
qchem -nt $SLURM_NTASKS qchem_2461.in qchem_2461.out
qchem -nt $SLURM_NTASKS qchem_2462.in qchem_2462.out
qchem -nt $SLURM_NTASKS qchem_2463.in qchem_2463.out
qchem -nt $SLURM_NTASKS qchem_2464.in qchem_2464.out
qchem -nt $SLURM_NTASKS qchem_2465.in qchem_2465.out
qchem -nt $SLURM_NTASKS qchem_2466.in qchem_2466.out
qchem -nt $SLURM_NTASKS qchem_2467.in qchem_2467.out
qchem -nt $SLURM_NTASKS qchem_2468.in qchem_2468.out
qchem -nt $SLURM_NTASKS qchem_2469.in qchem_2469.out
qchem -nt $SLURM_NTASKS qchem_2470.in qchem_2470.out
qchem -nt $SLURM_NTASKS qchem_2471.in qchem_2471.out
qchem -nt $SLURM_NTASKS qchem_2472.in qchem_2472.out
qchem -nt $SLURM_NTASKS qchem_2473.in qchem_2473.out
qchem -nt $SLURM_NTASKS qchem_2474.in qchem_2474.out
qchem -nt $SLURM_NTASKS qchem_2475.in qchem_2475.out
qchem -nt $SLURM_NTASKS qchem_2476.in qchem_2476.out
qchem -nt $SLURM_NTASKS qchem_2477.in qchem_2477.out
qchem -nt $SLURM_NTASKS qchem_2478.in qchem_2478.out
qchem -nt $SLURM_NTASKS qchem_2479.in qchem_2479.out
qchem -nt $SLURM_NTASKS qchem_2480.in qchem_2480.out
qchem -nt $SLURM_NTASKS qchem_2481.in qchem_2481.out
qchem -nt $SLURM_NTASKS qchem_2482.in qchem_2482.out
qchem -nt $SLURM_NTASKS qchem_2483.in qchem_2483.out
qchem -nt $SLURM_NTASKS qchem_2484.in qchem_2484.out
qchem -nt $SLURM_NTASKS qchem_2485.in qchem_2485.out
qchem -nt $SLURM_NTASKS qchem_2486.in qchem_2486.out
qchem -nt $SLURM_NTASKS qchem_2487.in qchem_2487.out
qchem -nt $SLURM_NTASKS qchem_2488.in qchem_2488.out
qchem -nt $SLURM_NTASKS qchem_2489.in qchem_2489.out
qchem -nt $SLURM_NTASKS qchem_2490.in qchem_2490.out
qchem -nt $SLURM_NTASKS qchem_2491.in qchem_2491.out
qchem -nt $SLURM_NTASKS qchem_2492.in qchem_2492.out
qchem -nt $SLURM_NTASKS qchem_2493.in qchem_2493.out
qchem -nt $SLURM_NTASKS qchem_2494.in qchem_2494.out
qchem -nt $SLURM_NTASKS qchem_2495.in qchem_2495.out
qchem -nt $SLURM_NTASKS qchem_2496.in qchem_2496.out
qchem -nt $SLURM_NTASKS qchem_2497.in qchem_2497.out
qchem -nt $SLURM_NTASKS qchem_2498.in qchem_2498.out
qchem -nt $SLURM_NTASKS qchem_2499.in qchem_2499.out
qchem -nt $SLURM_NTASKS qchem_2500.in qchem_2500.out
qchem -nt $SLURM_NTASKS qchem_2501.in qchem_2501.out
qchem -nt $SLURM_NTASKS qchem_2502.in qchem_2502.out
qchem -nt $SLURM_NTASKS qchem_2503.in qchem_2503.out
qchem -nt $SLURM_NTASKS qchem_2504.in qchem_2504.out
qchem -nt $SLURM_NTASKS qchem_2505.in qchem_2505.out
qchem -nt $SLURM_NTASKS qchem_2506.in qchem_2506.out
qchem -nt $SLURM_NTASKS qchem_2507.in qchem_2507.out
qchem -nt $SLURM_NTASKS qchem_2508.in qchem_2508.out
qchem -nt $SLURM_NTASKS qchem_2509.in qchem_2509.out
qchem -nt $SLURM_NTASKS qchem_2510.in qchem_2510.out
qchem -nt $SLURM_NTASKS qchem_2511.in qchem_2511.out
qchem -nt $SLURM_NTASKS qchem_2512.in qchem_2512.out
qchem -nt $SLURM_NTASKS qchem_2513.in qchem_2513.out
qchem -nt $SLURM_NTASKS qchem_2514.in qchem_2514.out
qchem -nt $SLURM_NTASKS qchem_2515.in qchem_2515.out
qchem -nt $SLURM_NTASKS qchem_2516.in qchem_2516.out
qchem -nt $SLURM_NTASKS qchem_2517.in qchem_2517.out
qchem -nt $SLURM_NTASKS qchem_2518.in qchem_2518.out
qchem -nt $SLURM_NTASKS qchem_2519.in qchem_2519.out
qchem -nt $SLURM_NTASKS qchem_2520.in qchem_2520.out
qchem -nt $SLURM_NTASKS qchem_2521.in qchem_2521.out
qchem -nt $SLURM_NTASKS qchem_2522.in qchem_2522.out
qchem -nt $SLURM_NTASKS qchem_2523.in qchem_2523.out
qchem -nt $SLURM_NTASKS qchem_2524.in qchem_2524.out
qchem -nt $SLURM_NTASKS qchem_2525.in qchem_2525.out
qchem -nt $SLURM_NTASKS qchem_2526.in qchem_2526.out
qchem -nt $SLURM_NTASKS qchem_2527.in qchem_2527.out
qchem -nt $SLURM_NTASKS qchem_2528.in qchem_2528.out
qchem -nt $SLURM_NTASKS qchem_2529.in qchem_2529.out
qchem -nt $SLURM_NTASKS qchem_2530.in qchem_2530.out
qchem -nt $SLURM_NTASKS qchem_2531.in qchem_2531.out
qchem -nt $SLURM_NTASKS qchem_2532.in qchem_2532.out
qchem -nt $SLURM_NTASKS qchem_2533.in qchem_2533.out
qchem -nt $SLURM_NTASKS qchem_2534.in qchem_2534.out
qchem -nt $SLURM_NTASKS qchem_2535.in qchem_2535.out
qchem -nt $SLURM_NTASKS qchem_2536.in qchem_2536.out
qchem -nt $SLURM_NTASKS qchem_2537.in qchem_2537.out
qchem -nt $SLURM_NTASKS qchem_2538.in qchem_2538.out
qchem -nt $SLURM_NTASKS qchem_2539.in qchem_2539.out
qchem -nt $SLURM_NTASKS qchem_2540.in qchem_2540.out
qchem -nt $SLURM_NTASKS qchem_2541.in qchem_2541.out
qchem -nt $SLURM_NTASKS qchem_2542.in qchem_2542.out
qchem -nt $SLURM_NTASKS qchem_2543.in qchem_2543.out
qchem -nt $SLURM_NTASKS qchem_2544.in qchem_2544.out
qchem -nt $SLURM_NTASKS qchem_2545.in qchem_2545.out
qchem -nt $SLURM_NTASKS qchem_2546.in qchem_2546.out
qchem -nt $SLURM_NTASKS qchem_2547.in qchem_2547.out
qchem -nt $SLURM_NTASKS qchem_2548.in qchem_2548.out
qchem -nt $SLURM_NTASKS qchem_2549.in qchem_2549.out
qchem -nt $SLURM_NTASKS qchem_2550.in qchem_2550.out
qchem -nt $SLURM_NTASKS qchem_2551.in qchem_2551.out
qchem -nt $SLURM_NTASKS qchem_2552.in qchem_2552.out
qchem -nt $SLURM_NTASKS qchem_2553.in qchem_2553.out
qchem -nt $SLURM_NTASKS qchem_2554.in qchem_2554.out
qchem -nt $SLURM_NTASKS qchem_2555.in qchem_2555.out
qchem -nt $SLURM_NTASKS qchem_2556.in qchem_2556.out
qchem -nt $SLURM_NTASKS qchem_2557.in qchem_2557.out
qchem -nt $SLURM_NTASKS qchem_2558.in qchem_2558.out
qchem -nt $SLURM_NTASKS qchem_2559.in qchem_2559.out
qchem -nt $SLURM_NTASKS qchem_2560.in qchem_2560.out
qchem -nt $SLURM_NTASKS qchem_2561.in qchem_2561.out
qchem -nt $SLURM_NTASKS qchem_2562.in qchem_2562.out
qchem -nt $SLURM_NTASKS qchem_2563.in qchem_2563.out
qchem -nt $SLURM_NTASKS qchem_2564.in qchem_2564.out
qchem -nt $SLURM_NTASKS qchem_2565.in qchem_2565.out
qchem -nt $SLURM_NTASKS qchem_2566.in qchem_2566.out
qchem -nt $SLURM_NTASKS qchem_2567.in qchem_2567.out
qchem -nt $SLURM_NTASKS qchem_2568.in qchem_2568.out
qchem -nt $SLURM_NTASKS qchem_2569.in qchem_2569.out
qchem -nt $SLURM_NTASKS qchem_2570.in qchem_2570.out
qchem -nt $SLURM_NTASKS qchem_2571.in qchem_2571.out
qchem -nt $SLURM_NTASKS qchem_2572.in qchem_2572.out
qchem -nt $SLURM_NTASKS qchem_2573.in qchem_2573.out
qchem -nt $SLURM_NTASKS qchem_2574.in qchem_2574.out
qchem -nt $SLURM_NTASKS qchem_2575.in qchem_2575.out
qchem -nt $SLURM_NTASKS qchem_2576.in qchem_2576.out
qchem -nt $SLURM_NTASKS qchem_2577.in qchem_2577.out
qchem -nt $SLURM_NTASKS qchem_2578.in qchem_2578.out
qchem -nt $SLURM_NTASKS qchem_2579.in qchem_2579.out
qchem -nt $SLURM_NTASKS qchem_2580.in qchem_2580.out
qchem -nt $SLURM_NTASKS qchem_2581.in qchem_2581.out
qchem -nt $SLURM_NTASKS qchem_2582.in qchem_2582.out
qchem -nt $SLURM_NTASKS qchem_2583.in qchem_2583.out
qchem -nt $SLURM_NTASKS qchem_2584.in qchem_2584.out
qchem -nt $SLURM_NTASKS qchem_2585.in qchem_2585.out
qchem -nt $SLURM_NTASKS qchem_2586.in qchem_2586.out
qchem -nt $SLURM_NTASKS qchem_2587.in qchem_2587.out
qchem -nt $SLURM_NTASKS qchem_2588.in qchem_2588.out
qchem -nt $SLURM_NTASKS qchem_2589.in qchem_2589.out
qchem -nt $SLURM_NTASKS qchem_2590.in qchem_2590.out
qchem -nt $SLURM_NTASKS qchem_2591.in qchem_2591.out
qchem -nt $SLURM_NTASKS qchem_2592.in qchem_2592.out
qchem -nt $SLURM_NTASKS qchem_2593.in qchem_2593.out
qchem -nt $SLURM_NTASKS qchem_2594.in qchem_2594.out
qchem -nt $SLURM_NTASKS qchem_2595.in qchem_2595.out
qchem -nt $SLURM_NTASKS qchem_2596.in qchem_2596.out
qchem -nt $SLURM_NTASKS qchem_2597.in qchem_2597.out
qchem -nt $SLURM_NTASKS qchem_2598.in qchem_2598.out
qchem -nt $SLURM_NTASKS qchem_2599.in qchem_2599.out
qchem -nt $SLURM_NTASKS qchem_2600.in qchem_2600.out
qchem -nt $SLURM_NTASKS qchem_2601.in qchem_2601.out
qchem -nt $SLURM_NTASKS qchem_2602.in qchem_2602.out
qchem -nt $SLURM_NTASKS qchem_2603.in qchem_2603.out
qchem -nt $SLURM_NTASKS qchem_2604.in qchem_2604.out
qchem -nt $SLURM_NTASKS qchem_2605.in qchem_2605.out
qchem -nt $SLURM_NTASKS qchem_2606.in qchem_2606.out
qchem -nt $SLURM_NTASKS qchem_2607.in qchem_2607.out
qchem -nt $SLURM_NTASKS qchem_2608.in qchem_2608.out
qchem -nt $SLURM_NTASKS qchem_2609.in qchem_2609.out
qchem -nt $SLURM_NTASKS qchem_2610.in qchem_2610.out
qchem -nt $SLURM_NTASKS qchem_2611.in qchem_2611.out
qchem -nt $SLURM_NTASKS qchem_2612.in qchem_2612.out
qchem -nt $SLURM_NTASKS qchem_2613.in qchem_2613.out
qchem -nt $SLURM_NTASKS qchem_2614.in qchem_2614.out
qchem -nt $SLURM_NTASKS qchem_2615.in qchem_2615.out
qchem -nt $SLURM_NTASKS qchem_2616.in qchem_2616.out
qchem -nt $SLURM_NTASKS qchem_2617.in qchem_2617.out
qchem -nt $SLURM_NTASKS qchem_2618.in qchem_2618.out
qchem -nt $SLURM_NTASKS qchem_2619.in qchem_2619.out
qchem -nt $SLURM_NTASKS qchem_2620.in qchem_2620.out
qchem -nt $SLURM_NTASKS qchem_2621.in qchem_2621.out
qchem -nt $SLURM_NTASKS qchem_2622.in qchem_2622.out
qchem -nt $SLURM_NTASKS qchem_2623.in qchem_2623.out
qchem -nt $SLURM_NTASKS qchem_2624.in qchem_2624.out
qchem -nt $SLURM_NTASKS qchem_2625.in qchem_2625.out
qchem -nt $SLURM_NTASKS qchem_2626.in qchem_2626.out
qchem -nt $SLURM_NTASKS qchem_2627.in qchem_2627.out
qchem -nt $SLURM_NTASKS qchem_2628.in qchem_2628.out
qchem -nt $SLURM_NTASKS qchem_2629.in qchem_2629.out
qchem -nt $SLURM_NTASKS qchem_2630.in qchem_2630.out
qchem -nt $SLURM_NTASKS qchem_2631.in qchem_2631.out
qchem -nt $SLURM_NTASKS qchem_2632.in qchem_2632.out
qchem -nt $SLURM_NTASKS qchem_2633.in qchem_2633.out
qchem -nt $SLURM_NTASKS qchem_2634.in qchem_2634.out
qchem -nt $SLURM_NTASKS qchem_2635.in qchem_2635.out
qchem -nt $SLURM_NTASKS qchem_2636.in qchem_2636.out
qchem -nt $SLURM_NTASKS qchem_2637.in qchem_2637.out
qchem -nt $SLURM_NTASKS qchem_2638.in qchem_2638.out
qchem -nt $SLURM_NTASKS qchem_2639.in qchem_2639.out
qchem -nt $SLURM_NTASKS qchem_2640.in qchem_2640.out
qchem -nt $SLURM_NTASKS qchem_2641.in qchem_2641.out
qchem -nt $SLURM_NTASKS qchem_2642.in qchem_2642.out
qchem -nt $SLURM_NTASKS qchem_2643.in qchem_2643.out
qchem -nt $SLURM_NTASKS qchem_2644.in qchem_2644.out
qchem -nt $SLURM_NTASKS qchem_2645.in qchem_2645.out
 
echo 'Finishing program on'
date
 
