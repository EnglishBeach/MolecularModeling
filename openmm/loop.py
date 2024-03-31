from pathlib import Path
import os
import simulator

for path in os.listdir('systems'):
    sim_path = Path(f'systems/{path}')
    substance, x, rho = path.split('_')

    print('#' * 30, f'{substance:10}: {int(x):<3}')
    for T in [15, 20, 25, 30, 40]:
        result_dir = Path(f'results') / sim_path.stem / str(T)
        if result_dir.exists():
            continue
        result_dir = result_dir.as_posix()

        print('#' * 20, 'Tempetature: ', T)
        sim = simulator.Simulation(T=T)
        sim.load(sim_path)
        sim.equilibrate(result_dir, steps=20)
        sim.start_product(result_dir)
        sim.simulate(500)
        sim.get_msd().to_csv(f'{result_dir}/msd.csv')
