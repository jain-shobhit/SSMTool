function [data y] = coco_funcfile(prob, data, u)
% The function coco_funcfile demonstrates the required
% function definition syntax for encoding a coco-compatible
% zero or monitor function. 
%
% Description of arguments:
% 1. The input argument prob contains a continuation problem structure
% initialized using the coco_prob utility.
% 2. The input and output argument data contains the function data
% structure whose fields may be used and modified in the function body.
% The initial content of this argument is assigned in the call to the
% coco_add_func constructor. 
% 3. The input argument u contains a vector of numerical values for
% the subset of continuation variables associated with the functions
% dependency index set, assigned in the call to the coco_add_func
% constructor.
% 4. The return argument y contains a vector of numerical values of
% the component zero or monitor functions.

% Example usage 1:
% Extract individual components of the input argument u and assign a
% functional expression to the return argument y. This approach encodes a
% specific functional form for the zero or monitor functions in terms of
% the components of u.
% y = [u(1)^2+(u(2)-1)^2-1; u(3)+2*u(2)-u(1)-3];

% Example usage 2:
% Extract individual components of the input argument u and assign a
% functional expression to the return argument y. This approach encodes a
% specific functional form for the zero or monitor functions in terms of
% the components of u and parameterized by a tolerance stored in the
% function data structure.
% fhan = @(x) sin(x)./x;
% y = quad(fhan, 0, u(1), data.tol) - u(2);

% Example usage 3:
% Use index sets stored in the function data structure to extract
% subsets of the input argument u. Assign the output of a call to a
% function whose handle is stored in the function data structure to the
% return argument y. This approach encodes a general functional form for
% the zero or monitor functions in terms of two distinct input arguments
% for a non-coco-compatible encoding of the zero or monitor functions.
% x = u(data.xidx);
% p = u(data.pidx);
% y = data.fhan(x, p);

% Example usage 4:
% Use index sets stored in the function data structure to extract
% subsets of the input argument u. Assign the outputs of a call to a
% function whose handle is stored in the function data structure to the
% return arguments data and y. This approach encodes a general functional
% form for the zero or monitor functions in terms of two distinct input
% arguments for a non-coco-compatible encoding of the zero or monitor
% functions, allowing for the possibility of changes to the function data
% structure.
% x = u(data.xidx);
% p = u(data.pidx);
% [data.opts y] = data.fhan(x, p, data.opts);

end