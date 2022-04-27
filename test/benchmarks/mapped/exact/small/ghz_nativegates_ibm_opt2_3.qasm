// i 0 1 2
// o 0 1 2
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
rz(1.5707963267949) q[2];
sx q[2];
rz(1.5707963267949) q[2];
cx q[2], q[1];
cx q[1], q[0];

