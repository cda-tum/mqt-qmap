OPENQASM 2.0;
include "qelib1.inc";
gate ccz q0,q1,q2 { p(pi/4) q1; p(pi/4) q2; cx q1,q2; u(0,0,-pi/4) q2; cx q1,q2; u(0,0,0) q2; cx q1,q0; p(-pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,pi/4) q2; cx q0,q2; u(0,-pi/2,0) q2; cx q1,q0; p(pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,-pi/4) q2; cx q0,q2; u(0,0,0) q2; }
gate ccz_139716245561296 q0,q1,q2 { p(pi/4) q1; p(pi/4) q2; cx q1,q2; u(0,0,-pi/4) q2; cx q1,q2; u(0,0,0) q2; cx q1,q0; p(-pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,pi/4) q2; cx q0,q2; u(0,-pi/2,0) q2; cx q1,q0; p(pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,-pi/4) q2; cx q0,q2; u(0,0,0) q2; }
gate ccz_o0 q0,q1,q2 { x q0; x q1; ccz_139716245561296 q0,q1,q2; x q0; x q1; }
gate ccz_139716245366384 q0,q1,q2 { p(pi/4) q1; p(pi/4) q2; cx q1,q2; u(0,0,-pi/4) q2; cx q1,q2; u(0,0,0) q2; cx q1,q0; p(-pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,pi/4) q2; cx q0,q2; u(0,-pi/2,0) q2; cx q1,q0; p(pi/4) q0; p(pi/4) q2; cx q0,q2; u(0,0,-pi/4) q2; cx q0,q2; u(0,0,0) q2; }
gate cz_o0 q0,q1 { x q0; cz q0,q1; x q0; }
gate mcx q0,q1,q2,q3 { h q3; p(pi/8) q0; p(pi/8) q1; p(pi/8) q2; p(pi/8) q3; cx q0,q1; p(-pi/8) q1; cx q0,q1; cx q1,q2; p(-pi/8) q2; cx q0,q2; p(pi/8) q2; cx q1,q2; p(-pi/8) q2; cx q0,q2; cx q2,q3; p(-pi/8) q3; cx q1,q3; p(pi/8) q3; cx q2,q3; p(-pi/8) q3; cx q0,q3; p(pi/8) q3; cx q2,q3; p(-pi/8) q3; cx q1,q3; p(pi/8) q3; cx q2,q3; p(-pi/8) q3; cx q0,q3; h q3; }
gate gate_Q q0,q1,q2,q3 { ccz q0,q2,q3; ccz_o0 q0,q3,q1; x q2; z q2; x q2; ccz_139716245366384 q0,q1,q2; cz_o0 q3,q2; h q3; h q2; h q1; h q0; x q0; x q1; x q2; x q3; h q3; mcx q0,q1,q2,q3; h q3; x q0; x q1; x q2; x q3; h q0; h q1; h q2; h q3; }
qreg q[4];
h q[0];
h q[1];
h q[2];
h q[3];
gate_Q q[0],q[1],q[2],q[3];
