%%
length = 21;
N = 1;
ang = [60];

for i = 1 : N
    PsfEstim{i} = motion_filter(length, ang(i));
end
save(['PsfEstim.mat'],'PsfEstim');%����ģ���˵ľ���

%% ����ģ���ˣ�
C = [];D = ones(length,3);
for k = 1:N 
    C = [C,D,PsfEstim{k}/max(PsfEstim{k}(:))]; 
end
figure, imshow(C , []);  title('all PSF normalized ')
imwrite(C,['/kernelsEstim.png']);