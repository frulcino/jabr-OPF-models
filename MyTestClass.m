classdef MyTestClass < matlab.unittest.TestCase
    methods (Test)
        function testSignFunction(testCase)
            mpc = loadcase("case9.m");
            result = [ss(1,4,mpc), ss(5,4,mpc),ss(9,2,mpc), ss(9,9,mpc)];
            expected = [1,-1,0,0];
            testCase.verifyEqual(result, expected);
        end
        
        function testIdside(testCase)
            mpc = loadcase("case9.m");
            result = [idside(1,4,mpc),idside(6,4,mpc),idside(9,4,mpc)];
            expected = [1,0,9];
            testCase.verifyEqual(result, expected);
        end
        
        function testNmat(testCase)
            mpc = loadcase("case9.m");
            result = Nmat(4,mpc);
            expected = sparse([1,2,9],[2,1,2],[1,1,1])
            testCase.verifyEqual(result, expected);
        end
    end
end