
function [index] = idside(a, b, mpc)
    % Given two connected buses a and b, return the branch index in mpc.branch
    size_branch = size(mpc.branch);
    Nbranch = size_branch(1);
    index = 0; % Initialize index to zero
    
    for i = 1:Nbranch
        f = mpc.branch(i, 1);
        t = mpc.branch(i, 2);
        
        if isequal([a, b], [f, t]) || isequal([b, a], [f, t])
            index = i;
            break; % Exit the loop once the index is found
        end
    end
end

