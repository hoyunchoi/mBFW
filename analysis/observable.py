from enum import Enum


class Observable(Enum):
    maximum_cluster_size = "maximum_cluster_size"
    second_maximum_cluster_size = "second_maximum_cluster_size"
    mean_cluster_size = "mean_cluster_size"
    second_moment = "second_moment"
    inter_event_time = "inter_event_time"
    max_delta_order_parameter = "max_delta_order_parameter"
