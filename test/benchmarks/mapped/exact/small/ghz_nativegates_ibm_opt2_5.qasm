// i 0 1 2 3 5
// o 0 1 2 3 5
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
rz(1.5707963267949) q[5];
sx q[5];
rz(1.5707963267949) q[5];
cx q[5], q[3];
cx q[3], q[2];
cx q[2], q[1];
cx q[1], q[0];

