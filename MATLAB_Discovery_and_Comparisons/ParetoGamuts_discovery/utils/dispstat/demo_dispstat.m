dispstat('','init'); %one time only init 

dispstat('Begining the process...','keepthis','timestamp'); 
for i = 1:100 
dispstat(sprintf('Processing %d%%',i),'timestamp'); 
dispstat('hi', 'keepprev', 'keepthis');
%doing some heavy stuff here 
pause(0.1)
end
dispstat('Finished.','keepprev');

dispstat('Begining the process...','keepthis','timespamp'); 
for i = 1:100 
dispstat(sprintf('Processing %d%%',i),'timestamp'); 
%doing some heavy stuff here 
pause(0.1)
end
dispstat('Finished.','keepprev');