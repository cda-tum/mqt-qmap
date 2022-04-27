// i 0 1
// o 0 1
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
x q[1];
cx q[1], q[0];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
cx q[0], q[1];

