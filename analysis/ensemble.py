import os
import numpy as np
import pathlib

from name import Name
from observable import Observable

class Ensemble:
    data_dir = "/pds/pds11/hoyun/mBFW/data/"

    def __init__(self, network_size: int, acceptance_threshold: float, remove: bool) -> None:
        self.network_size = network_size
        self.acceptance_threshold = acceptance_threshold
        self.remove = remove

        self.name = Name(self.network_size, self.acceptance_threshold)

    @staticmethod
    def __extract_ensemble_size_from_file_name(file_name: str) -> int:
        file_name_list = file_name.replace('-', ',').replace('.', ',').split(',')
        for fraction in file_name_list:
            if "E" in fraction:
                return int(fraction.removeprefix("E"))

        raise ValueError(f"No ensemble prefix at {file_name}")

    def __scan_directory(self, directory: str) -> tuple[list[str]]:
        # *Get list of file inside directory
        file_name_list = [file_name
                          for file_name in os.listdir(directory)
                          if self.name.NG in file_name]

        # *Get list of ensemble sizes from file name list
        ensemble_list = [self.__extract_ensemble_size_from_file_name(file_name)
                         for file_name in file_name_list]

        #* Get full path of file name
        file_path_list = [os.path.join(directory, file_name)
                          for file_name in file_name_list]
        return file_path_list, ensemble_list

    def __remove_files(self, file_path_list:list[str]) -> None:
        # Remove files at file_path_list if self.remove
        if self.remove:
            for file_path in file_path_list:
                os.remove(file_path)

    @staticmethod
    def __additional_directory(parent: str, directory: str) -> str:
        new_directory = os.path.join(parent, directory)
        pathlib.Path(new_directory).mkdir(parents=True, exist_ok=True)
        return new_directory

    def plain(self, observable: Observable) -> None:
        directory = os.path.join(Ensemble.data_dir, observable.value)
        file_path_list, ensemble_list = self.__scan_directory(directory)
        if len(file_path_list) == 1:
            return
        tot_ensemble = sum(ensemble_list)

        avg_data = np.concatenate([np.loadtxt(file_path)[np.newaxis, :] for file_path in file_path_list])
        avg_data = np.average(avg_data, axis=0, weights=ensemble_list)

        #* Remove previous data
        self.__remove_files(file_path_list)

        #* Save plain text
        np.savetxt(os.path.join(directory, self.name.get_NGE(tot_ensemble)), avg_data,
                   delimiter=",",
                   fmt="%.10f")

        #* Save npy data
        ensemble_dir = self.__additional_directory(directory, "ensemble")
        np.save(os.path.join(ensemble_dir, self.name.get_NGE(tot_ensemble, npy=True)),
                avg_data)

    def inter_event_time(self) -> None:
        directory = os.path.join(Ensemble.data_dir, Observable.inter_event_time.value)
        file_path_list, ensemble_list = self.__scan_directory(directory)
        if len(file_path_list) == 1:
            return
        tot_ensemble = sum(ensemble_list)

        avg_data = dict()
        for file_path in file_path_list:
            data = np.loadtxt(file_path, delimiter=",")
            for single in data:
                time = single[0]
                try:
                    avg_data[time][0] = np.average([avg_data[time][0], single[1]],
                                                   weights=[avg_data[time][1], single[2]])
                    avg_data[time][1] = avg_data[time][1] + single[2]
                except KeyError:
                    avg_data[time] = single[1:]
        avg_data = np.array([np.insert(val, 0, key) for (key,val) in avg_data.items()])

        #* Remove previous data
        self.__remove_files(file_path_list)

        #* Save plain text
        np.savetxt(os.path.join(directory, self.name.get_NGE(tot_ensemble)), avg_data,
                   delimiter=",",
                   fmt="%.10f")

        #* Save npy data
        ensemble_dir = self.__additional_directory(directory, "ensemble")
        np.save(os.path.join(ensemble_dir, self.name.get_NGE(tot_ensemble, npy=True)),
                avg_data[:, :2]/self.network_size)

    def max_delta_order_parameter(self) -> None:
        directory = os.path.join(Ensemble.data_dir, Observable.max_delta_order_parameter.value)
        file_path_list, ensemble_list = self.__scan_directory(directory)
        if len(file_path_list) == 1:
            return
        tot_ensemble = sum(ensemble_list)

        tot_data = np.concatenate([np.loadtxt(file_path, delimiter=',') for file_path in file_path_list])

        #* Remove previous data
        self.__remove_files(file_path_list)

        #* Save plain text
        np.savetxt(os.path.join(directory, self.name.get_NGE(tot_ensemble)), tot_data,
                   delimiter=",",
                   fmt="%.10f")

        #* Save npy data
        ensemble_dir = self.__additional_directory(directory, "ensemble")
        np.save(os.path.join(ensemble_dir, self.name.get_NGE(tot_ensemble, npy=True)),
                tot_data)


def main():
    ensemble = Ensemble(10000, 0.5, remove=True)
    ensemble.plain(Observable.maximum_cluster_size)
    ensemble.plain(Observable.second_maximum_cluster_size)
    ensemble.plain(Observable.mean_cluster_size)
    ensemble.plain(Observable.second_moment)

    ensemble.inter_event_time()
    ensemble.max_delta_order_parameter()


if __name__ == "__main__":
    main()
