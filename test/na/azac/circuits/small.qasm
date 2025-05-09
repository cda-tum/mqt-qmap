OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
y q;
cz q[0],q[1];
cz q[0],q[2];
z q[1];
z q[2];
