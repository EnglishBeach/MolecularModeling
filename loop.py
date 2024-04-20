import time
from pathlib import Path

import base
from gromacs import simulator

for substance, combinations in base.DATAS.items():
    for x, rho in combinations:
        T = 25
        title = f'{substance.name :10}: {int(x):<3} - {T} C'
        workdir = Path(f'solvs_long/{substance.name}_{x}_{rho}')
        print(
            f"{title:*^100}",
        )
        if workdir.exists():
            continue
        sim = simulator.Simulation(
            rho=rho / 1000,
            T=T,
            x=x / 100,
            substance=substance.name,
            workdir=workdir,
        )
        sim.create_box()
        sim.center()

        time.sleep(2)
        sim.em()

        time.sleep(2)
        sim.nvt()

        time.sleep(2)
        sim.npt()

        time.sleep(2)
        sim.md()

        time.sleep(2)
        sim.summarize()
