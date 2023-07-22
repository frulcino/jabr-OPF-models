function [N] = Nmat(b,mpc)
%computes the neighbour matrix of b
% return a 2*Nbranch matrix N where N(1,l)=1 (N(2,l)=1) if b is the starting (arriving) bus of l
% 0 otherwise and 
   size_branch = size(mpc.branch);
   Nbranch = size_branch(1);
   N = sparse([],[],[],Nbranch,2);
   for i = 1:Nbranch
       f = mpc.branch(i,1);
       t = mpc.branch(i,2);
       if isequal(b,f);
           N(i,1) = 1;
       elseif isequal(b,t);
           N(i,2) = 1;
       end
   end
end

