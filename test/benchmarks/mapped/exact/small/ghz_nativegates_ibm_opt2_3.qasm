// i 2 1 0
// o 2 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[3];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
cx q[0], q[1];
cx q[1], q[2];
measure q[0] -> c[2];
measure q[1] -> c[1];
measure q[2] -> c[0];

