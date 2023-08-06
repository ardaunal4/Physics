ma = 1/2;
mb = 1/4;

w = sym('w');

Ma = [-w^2*ma+2 -1
         1 0];
Mb = [-w^2*mb+2 -1
         1 0];

ABA = Ma*Mb*Ma;
BAABA = Mb*Ma*Ma*Mb*Ma;
ABABAABA = Ma*Mb*Ma*Mb*Ma*Ma*Mb*Ma;
BAABAABABAABA = Mb*Ma*Ma*Mb*Ma*Ma*Mb*Ma*Mb*Ma*Ma*Mb*Ma;

trABA = trace(ABA);
equation1 = trABA == 2;
roots1 = vpasolve(equation1, w);
roots1 = simplify(roots1);

trBAABA = trace(BAABA);
equation2 = trBAABA == 2;
roots2 = vpasolve(equation2, w);
roots2 = simplify(roots2);

trABABAABA = trace(ABABAABA);
equation3 = trABABAABA == 2;
roots3 = vpasolve(equation3, w);
roots3 = simplify(roots3);

trBAABAABABAABA = trace(BAABAABABAABA);
equation4 = trBAABAABABAABA == 2;
roots4 = vpasolve(equation4, w);
roots4 = simplify(roots4);

