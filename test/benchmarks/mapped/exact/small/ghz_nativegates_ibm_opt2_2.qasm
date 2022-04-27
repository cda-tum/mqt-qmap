// i 0 1
// o 0 1
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
rz(1.5707963267949) q[1];
sx q[1];
rz(1.5707963267949) q[1];
cx q[1], q[0];

