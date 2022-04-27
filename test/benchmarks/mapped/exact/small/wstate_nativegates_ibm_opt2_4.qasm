// i 3 2 1 0
// o 3 2 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
sx q[3];
rz(0.785398163397448) q[3];
sx q[3];
sx q[2];
rz(0.61547971) q[2];
sx q[2];
sx q[1];
rz(0.523598775598299) q[1];
sx q[1];
x q[0];
cx q[0], q[1];
sx q[1];
rz(0.523598775598299) q[1];
sx q[1];
cx q[1], q[2];
sx q[2];
rz(0.61547971) q[2];
sx q[2];
cx q[2], q[3];
sx q[3];
rz(0.785398163397448) q[3];
sx q[3];
cx q[1], q[0];
cx q[2], q[1];
cx q[3], q[2];

