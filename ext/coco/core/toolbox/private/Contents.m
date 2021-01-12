
%Private functions of COCO.
%
%COCO options structure:
%   coco_set_defaults     - Initialise COCO options for continuation.
%
%Finite state machine (continuation method):










%
%Tangent predictor:



%
%Extended system:


%
%Bifurcation diagram:




%
%Saving data:


%
%Screen output:


%   nwtn_print_data       - Print information about Newton's iteration.
%   nwtn_print_headline   - Print headline for Newton's iteration.
%
%Default functions and event handlers:
%   coco_default_callback - Default handler for coco's call-back messages.
%   coco_default_update   - Default handler for coco's update message.
%   coco_default_linsolve - Default implementation for solving Ax=b.
%   coco_default_print    - Default implementation for printing additional data.
%   coco_default_save     - Default implementation listing additional save data.
%   efunc_DFDX            - Evaluate Jacobian of extended system at UP=[U;P].
%   efunc_F               - Evaluate extended system at UP=[U;P].
%   efunc_events_F        - Evaluate event function..
%   efunc_init            - Initialise extended system.
%   efunc_monitor_F       - Evaluate monitor functions.
