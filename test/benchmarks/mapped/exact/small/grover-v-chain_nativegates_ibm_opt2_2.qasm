// i 0 1
// o 0 1
OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
qreg q[27];
creg meas[2];
rz(1.5707963267949) q[0];
sx q[0];
rz(1.5707963267949) q[0];
x flag[0];
barrier q[0];
barrier flag[0];
measure q -> meas;
measure flag -> meas;

