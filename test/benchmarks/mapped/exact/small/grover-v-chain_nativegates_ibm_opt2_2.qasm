// i 0 1
// o 0 1
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg meas[2];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
x q[1];
barrier q;
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[0] -> meas[0];
measure q[1] -> meas[1];

