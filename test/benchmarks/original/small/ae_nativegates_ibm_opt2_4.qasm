// Benchmark was created by MQT Bench on 2022-04-07
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.1.0
// Qiskit version: {'qiskit-terra': '0.19.2', 'qiskit-aer': '0.10.3', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.18.3', 'qiskit-aqua': '0.9.5', 'qiskit': '0.34.2', 'qiskit-nature': '0.3.1', 'qiskit-finance': '0.3.1', 'qiskit-optimization': '0.3.1', 'qiskit-machine-learning': '0.3.1'}
// Used Gate Set: ['id', 'rz', 'sx', 'x', 'cx', 'reset']
// Optimization Level: 2

OPENQASM 2.0;
include "qelib1.inc";
qreg eval[3];
qreg q[1];
creg meas[4];
rz(pi/2) eval[0];
sx eval[0];
rz(5*pi/8) eval[0];
rz(pi/2) eval[1];
sx eval[1];
rz(3*pi/4) eval[1];
rz(pi/2) eval[2];
sx eval[2];
rz(pi) eval[2];
rz(-pi) q[0];
sx q[0];
rz(2.2142974) q[0];
sx q[0];
cx eval[0],q[0];
sx q[0];
rz(2.2142974) q[0];
sx q[0];
rz(-pi) q[0];
cx eval[0],q[0];
rz(-pi) q[0];
sx q[0];
rz(2.2142974) q[0];
sx q[0];
cx eval[1],q[0];
sx q[0];
rz(1.2870022) q[0];
sx q[0];
rz(-pi) q[0];
cx eval[1],q[0];
rz(-pi) q[0];
sx q[0];
rz(1.2870022) q[0];
sx q[0];
cx eval[2],q[0];
rz(-pi) q[0];
sx q[0];
rz(0.56758822) q[0];
sx q[0];
cx eval[2],q[0];
sx q[0];
rz(0.56758822) q[0];
sx q[0];
rz(-pi) q[0];
sx eval[2];
rz(pi/2) eval[2];
cx eval[1],eval[2];
rz(pi/4) eval[2];
cx eval[1],eval[2];
sx eval[1];
rz(pi/2) eval[1];
rz(-pi/4) eval[2];
cx eval[0],eval[2];
rz(pi/8) eval[2];
cx eval[0],eval[2];
cx eval[0],eval[1];
rz(pi/4) eval[1];
cx eval[0],eval[1];
sx eval[0];
rz(pi/2) eval[0];
rz(-pi/4) eval[1];
rz(-pi/8) eval[2];
barrier eval[0],eval[1],eval[2],q[0];
measure eval[0] -> meas[0];
measure eval[1] -> meas[1];
measure eval[2] -> meas[2];
measure q[0] -> meas[3];
