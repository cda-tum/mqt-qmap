from qiskit import QuantumCircuit

class NAComputation:
    def __init__(self) -> None:
        """This class represents a neutral atom computation. The code can be retrieved as a string."""

class NAStatePreparationSolver:
    def __init__(self) -> None: ...

    # todo: add documentation
    # todo: verify that inputs are actually uint16
    def init(
        self,
        new_max_x: int,
        new_max_y: int,
        new_max_c: int,
        new_max_r: int,
        new_max_h_offset: int,
        new_max_v_offset: int,
        new_max_h_dist: int,
        new_max_vdist: int,
        new_min_entangling_y: int,
        new_max_entangling_y: int,
    ) -> None: ...

    class Result:
        def __init__(self) -> None: ...
        def yaml(self, indent: int = ..., compact: bool = ...) -> str: ...

    def solve(
        self,
        ops: list[tuple[int, int]],
        new_num_qubits: int,
        new_num_stages: int,
        new_num_transfers: int | None = ...,
        min_ops_order: bool = ...,
        shield_idle_qubits: bool = ...,
    ) -> Result: ...

def get_ops_for_solver(circ: QuantumCircuit) -> list[tuple[int, int]]: ...
def generate_code(
    circ: object,
    result: NAStatePreparationSolver.Result,
    max_hoffset: int,
    max_v_offset: int,
    min_entangling_y: int,
    max_entangling: int,
) -> NAComputation: ...
