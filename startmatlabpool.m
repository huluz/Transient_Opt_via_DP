function startmatlabpool(size)
	isstart = 0;
	nlabs = matlabpool('size');

	if nlabs == 0
		isstart = 1;
	end

	if isstart == 1
		if nargin == 0
			matlabpool('open','local');
		else
			try
				matlabpool('open','local',size);
			catch ce
				matlabpool('open','local');
				size = matlabpool('size');
				display(ce.message);
				display(strcat('输入size不正确，采用默认配置，size = ',num2str(size)));
			end
		end
	else
		display('matlabpool 已启动');
		if nargin == 0
			return;
			%matlabpool close;
			%startmatlabpool();
		elseif  nlabs ~= size
			matlabpool close;
			startmatlabpool(size);
		end
	end
