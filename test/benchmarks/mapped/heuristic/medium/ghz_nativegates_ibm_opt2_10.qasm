// i 13 14 11 9 8 5 3 2 1 0 10 7 12 6 4 15 16 17 18 19 20 21 22 23 24 25 26
// o 13 14 11 9 8 5 3 2 1 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[10];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
cx q[0], q[1];
cx q[1], q[2];
cx q[2], q[3];
cx q[3], q[5];
cx q[5], q[8];
cx q[8], q[9];
cx q[9], q[8];
cx q[8], q[11];
cx q[11], q[14];
cx q[14], q[13];
measure q[0] -> c[9];
measure q[1] -> c[8];
measure q[2] -> c[7];
measure q[3] -> c[6];
measure q[5] -> c[5];
measure q[8] -> c[4];
measure q[9] -> c[3];
measure q[11] -> c[2];
measure q[13] -> c[0];
measure q[14] -> c[1];
