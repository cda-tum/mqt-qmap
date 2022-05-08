// i 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
// o 0 1
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[2];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
x q[1];
cx q[1], q[0];
sx q[0];
rz(0.785398163397448) q[0];
sx q[0];
cx q[0], q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
