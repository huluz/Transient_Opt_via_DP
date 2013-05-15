%批量计算优化方案脚本

%上升情况
clc;clear;
load Increasing_InitialData;
TTOviaDP_Parallel_re03(ff,time,dt,qsmin,qsmax,dq,qbasic,fn);

%下降情况
load Decreasing_InitialData;
TTOviaDP_Parallel_re03(ff,time,dt,qsmin,qsmax,dq,qbasic,fn);

%波动情况
load Ramping_InitialData;
TTOviaDP_Parallel_re03(ff,time,dt,qsmin,qsmax,dq,qbasic,fn);

%长时间缓波动情况
load Slow_Ramping_24H_InitialData;
TTOviaDP_Parallel_re03(ff,time,dt,qsmin,qsmax,dq,qbasic,fn);

%长时间剧烈波动情况
load Severe_Ramping_24H_InitialData;
TTOviaDP_Parallel_re03(ff,time,dt,qsmin,qsmax,dq,qbasic,fn);