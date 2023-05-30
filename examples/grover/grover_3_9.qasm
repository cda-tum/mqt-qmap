OPENQASM 2.0;
include "qelib1.inc";
gate gate_Q q0,q1,q2 { cz q1,q2; x q1; z q1; x q1; h q2; h q1; h q0; x q0; x q1; x q2; h q2; ccx q0,q1,q2; h q2; x q0; x q1; x q2; h q0; h q1; h q2; }
qreg q[3];
h q[0];
h q[1];
h q[2];
gate_Q q[0],q[1],q[2];
