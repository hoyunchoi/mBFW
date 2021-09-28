class Name:
    def __init__(self, network_size: int, acceptance_threshold: float) -> None:
        self.network_size = network_size
        self.acceptance_threshold = acceptance_threshold

        self.NG = f'N{network_size},G{acceptance_threshold:.1f}'

    def get_NG(self, core_num: int=0, npy=False) -> str:
        if npy:
            return self.NG + ".npy"
        else:
            return self.NG + f"-{core_num}.dat"

    def get_NGE(self, ensemble_size: int, core_num: int=0, npy=False) -> str:
        if npy:
            return self.NG + f',E{ensemble_size}.npy'
        else:
            return self.NG + f',E{ensemble_size}-{core_num}.dat'
