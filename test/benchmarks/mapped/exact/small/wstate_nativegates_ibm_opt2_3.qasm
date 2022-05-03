// i 0 1 2
// o 0 1 2
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[3];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
sx q[1];
rz(0.61547971) q[1];
sx q[1];
x q[2];
cx q[2], q[1];
sx q[1];
rz(0.61547971) q[1];
sx q[1];
cx q[1], q[0];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
cx q[1], q[2];
cx q[0], q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];

