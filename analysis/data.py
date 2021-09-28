import os
import numpy as np

from name import Name
from observable import Observable


class Data:
    data_dir = "/pds/pds11/hoyun/mBFW/data/"

    def __init__(self, network_size: int, acceptance_threshold: float) -> None:
        self.network_size = network_size
        self.acceptance_threshold = acceptance_threshold
        self.name = Name(self.network_size, self.acceptance_threshold)

        self.time = np.linspace(0.0, 1.0, self.network_size)

    def read(self, observable_list: list[Observable]):
        if Observable.maximum_cluster_size in observable_list:
            self.maximum_cluster_size = self.__plain(Observable.maximum_cluster_size)

        if Observable.second_maximum_cluster_size in observable_list:
            self.second_maximum_cluster_size = self.__plain(Observable.second_maximum_cluster_size)

        if Observable.mean_cluster_size in observable_list:
            self.mean_cluster_size = self.__plain(Observable.mean_cluster_size)

        if Observable.second_moment in observable_list:
            self.second_moment = self.__plain(Observable.second_moment)

        if Observable.inter_event_time in observable_list:
            self.inter_event_time_t, self.inter_event_time_avg = self.__inter_event_time()

        if Observable.max_delta_order_parameter in observable_list:
            self.max_delta_order_parameter_t, self.max_delta_order_parameter_avg = self.__max_delta_order_parameter()

    def __get_ensemble_file_path(self, observable: Observable) -> str:
        directory = os.path.join(Data.data_dir, observable.value, "ensemble")
        file_name_list = [file_name
                          for file_name in os.listdir(directory)
                          if self.name.NG in file_name]
        if len(file_name_list) != 1:
            raise Exception(f"Check {directory}: N{self.network_size},G{self.acceptance_threshold}")

        return os.path.join(directory, file_name_list[0])

    def __plain(self, observable: str) -> np.ndarray:
        file_path = self.__get_ensemble_file_path(observable)
        return np.load(file_path)

    def __inter_event_time(self) -> tuple[np.ndarray]:
        file_path = self.__get_ensemble_file_path(Observable.inter_event_time)
        data = np.load(file_path)
        return data[:, 0], data[:, 1]

    def __max_delta_order_parameter(self) -> tuple[float]:
        file_path = self.__get_ensemble_file_path(Observable.inter_event_time)
        data = np.load(file_path)
        return np.average(data[:, 0]), np.average(data[:, 1])


if __name__ == "__main__":
    print("This is module data")
