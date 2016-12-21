

%%
varargin = {'a', 1, 'b', 2, 'pi', 3.14};

%% Standard way of operation:
clear P
P.a = [];
P.b = '0';
P.c = '3';
P.pi = pi;
P = util.parseInputs(P, varargin, 'error')

%% Should throw an error:
clear P
P.pi = 3;
P = util.parseInputs(P, varargin, 'error')

%% Add new values to the previously existing P:
clear P
P.c = '3';
P.pi = 3;
P
P = util.parseInputs(P, varargin, 'merge')

%% Should change values common to P and varargin, return others as R:
clear P
P.pi = 3;
[P, R] = util.parseInputs(P, varargin, 'split')


%% Give a structure as input:
clear Pin P
Pin.a = 1;
Pin.b = 2;
Pin.pi = 3.14;

P.a = 'a2';
P.c = 'c2';
P
P = util.parseInputs(P, Pin, 'merge')

%% Don't specify anything:
clear Pin P
Pin.a = 1;
Pin.b = 2;
Pin.pi = 3.14;
varargin = {Pin, 'c', 3}

P = util.parseInputs(struct(), varargin, 'merge')




