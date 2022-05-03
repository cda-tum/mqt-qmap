// i 3 2 1 0
// o 3 2 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[4];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
cx q[0], q[1];
cx q[1], q[2];
cx q[2], q[3];
measure q[0] -> c[3];
measure q[1] -> c[2];
measure q[2] -> c[1];
measure q[3] -> c[0];

