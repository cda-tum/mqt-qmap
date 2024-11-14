from typing import overload

from qiskit import QuantumCircuit

from mqt.core.na import NAComputation

Gate = Solver.Result.Gate
Qubit = Solver.Result.Qubit
Result = Solver.Result

class Solver:
    @overload
    def __init__(self) -> None: ...

    # todo: add documentation
    # todo: verify that inputs are actually uint16
    @overload
    def init(
        self,
        newMaxX: int,
        newMaxY: int,
        newMaxC: int,
        newMaxR: int,
        newMaxHOffset: int,
        newMaxVOffset: int,
        newMaxHDist: int,
        newMaxVdist: int,
        newMinEntanglingY: int,
        newMaxEntanglingY: int,
    ) -> None: ...

    class Result:
        class Qubit:
            @overload
            def __init__(self) -> None: ...
            @overload
            def __init__(self, x: int, y: int, a: bool, c: int, r: int, h: int, v: int) -> None: ...
            @overload
            def getX(self) -> int: ...
            @overload
            def getY(self) -> int: ...
            @overload
            def isAOD(self) -> bool: ...
            @overload
            def getC(self) -> int: ...
            @overload
            def getR(self) -> int: ...
            @overload
            def getH(self) -> int: ...
            @overload
            def getV(self) -> int: ...

        class Gate:
            @overload
            def __init__(self) -> None: ...
            @overload
            def __init__(self, stage: int, qubits: tuple[int, int]) -> None: ...
            @overload
            def getStage(self) -> int: ...
            @overload
            def getQubits(self) -> tuple[int, int]: ...

        class Stage:
            @overload
            def __init__(self) -> None: ...
            @overload
            def __init__(self, rydberg: bool, qubits: list[Qubit], gates: list[Gate]) -> None: ...
            @overload
            def isRydberg(self) -> bool: ...
            @overload
            def getQubits(self) -> list[Qubit]: ...
            @overload
            def getGates(self) -> list[Gate]: ...

        @overload
        def __init__(self) -> None: ...
        @overload
        def __init__(self, sat: bool) -> None: ...
        @overload
        def __init__(self, sat: bool, stages: list[Stage]) -> None: ...
        @overload
        def getStage(self, i: int) -> Stage: ...
        @overload
        def numStages(self) -> int: ...
        @overload
        def isSat(self) -> bool: ...
        @overload
        def yaml(self, indent: int = ..., compact: bool = ...) -> str: ...

    @overload
    def solve(
        self,
        ops: list[tuple[int, int]],
        newNumQubits: int,
        newNumStages: int,
        newNumTransfers: int = ...,
        minOpsOrder: bool = ...,
        shieldIdleQubits: bool = ...,
    ) -> Result: ...
    @overload
    def solve(
        self,
        ops: list[tuple[int, int]],
        newNumQubits: int,
        newNumStages: int,
        minOpsOrder: bool = ...,
        shieldIdleQubits: bool = ...,
    ) -> Result: ...

@overload
def getOpsForSolver(circ: QuantumCircuit) -> list[tuple[int, int]]: ...
@overload
def generate(
    circ: QuantumCircuit, result: Result, maxHoffset: int, maxVOffset: int, minEntanglingY: int, maxEntangling: int
) -> NAComputation: ...
