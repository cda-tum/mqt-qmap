// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// Qiskit version: {'qiskit-terra': '0.24.1', 'qiskit-aer': '0.12.0', 'qiskit-ignis': None, 'qiskit-ibmq-provider': '0.20.2', 'qiskit': '0.43.1', 'qiskit-nature': '0.6.2', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.5.0', 'qiskit-machine-learning': '0.6.1'}
// Used Gate Set: ['rx', 'rz', 'cz', 'measure']

OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[9];
rz(-pi/2) q[0];
rx(-pi/2) q[0];
rz(pi/2) q[1];
rx(-pi/2) q[1];
rz(pi/2) q[2];
rx(pi/2) q[2];
rz(pi/2) q[3];
rx(-pi/2) q[3];
rz(pi/2) q[4];
rx(pi/2) q[4];
rz(pi/2) q[5];
rx(-pi/2) q[5];
rz(pi/2) q[6];
rx(-pi/2) q[6];
rz(pi/2) q[7];
rx(pi/2) q[7];
rz(pi/2) q[8];
rx(-pi/2) q[8];
rx(pi) q[9];
cz q[0],q[9];
rx(pi/2) q[0];
rz(pi/2) q[0];
rz(-3*pi/2) q[9];
cz q[1],q[9];
rx(pi/2) q[1];
rz(-pi/2) q[1];
cz q[2],q[9];
rx(-pi/2) q[2];
rz(-pi/2) q[2];
cz q[3],q[9];
rx(pi/2) q[3];
rz(-pi/2) q[3];
cz q[4],q[9];
rx(-pi/2) q[4];
rz(-pi/2) q[4];
cz q[5],q[9];
rx(pi/2) q[5];
rz(-pi/2) q[5];
cz q[6],q[9];
rx(pi/2) q[6];
rz(-pi/2) q[6];
cz q[7],q[9];
rx(-pi/2) q[7];
rz(-pi/2) q[7];
cz q[8],q[9];
rx(pi/2) q[8];
rz(-pi/2) q[8];
rx(-pi/2) q[9];
rz(-pi/2) q[9];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
