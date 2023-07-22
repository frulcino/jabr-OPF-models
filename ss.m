function [sign] = bs(a,b,mpc)
% +1 is (a,b) is in branch, -1 if (b,a) is in branch and zero otherwise
   size_branch = size(mpc.branch);
   Nbranch = size_branch(1);
   sign = 0;
   for i = 1:Nbranch
       f = mpc.branch(i,1);
       t = mpc.branch(i,2);
       if isequal([a,b],[f,t]);
           sign = 1;
       elseif isequal([b,a],[f,t]);
           sign = -1;
       end
   end
end
