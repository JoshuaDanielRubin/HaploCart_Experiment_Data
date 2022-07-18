import bdsg
from bdsg.bdsg import ODGI
import sys

def f(x):
    graph.append_step(x, graph.get_handle_of_step(graph.path_begin(x)))
    graph.set_circularity(x, True)
    print(graph.get_is_circular(x))
    return True


graph = ODGI()
graph.deserialize(sys.argv[1])
graph.for_each_path_handle(lambda x: f(x))

graph.serialize(sys.argv[2])
