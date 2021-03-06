function [X]=generate_mo(name,opt)
% autogenerated by xsd_to_mo

if(nargin<2)
	opt=false;
end
version=1;
switch(version)
	case 1
		if(opt)
			X=generate_mo_withopt_v1(name);
		else
			X=generate_mo_noopt_v1(name);
		end
	otherwise
		error('Unsupported version number')
end
