OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
ccx q[2], q[4], q[6];
ccx q[2], q[5], q[7];
ccx q[3], q[4], q[7];
cx q[6], q[0];
cx q[7], q[1];
