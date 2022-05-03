// i 5 3 2 1 0 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
// o 5 3 2 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[5];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
cx q[0], q[1];
cx q[1], q[2];
cx q[2], q[3];
cx q[3], q[5];
measure q[0] -> c[4];
measure q[1] -> c[3];
measure q[2] -> c[2];
measure q[3] -> c[1];
measure q[5] -> c[0];
