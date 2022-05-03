// i 1 0
// o 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[2];
sx q[1];
rz(0.785398163397448) q[1];
sx q[1];
x q[0];
cx q[0], q[1];
sx q[1];
rz(0.785398163397448) q[1];
sx q[1];
cx q[1], q[0];
measure q[0] -> c[1];
measure q[1] -> c[0];

