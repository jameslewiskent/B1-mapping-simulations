function [ReorderPE1,ReorderPE2] = Calc_Reordering(settings)
% Calculate the reorder indices for both image trains. 
% We need to determine the ordering of the data, taking into account
% any partial Fourier or reduced phase resolution.

% First reorder PE1
% Determine which line is the centre (i.e. not the centre of Matrix_Size for partial Fourier)
Centre_Line = ceil(settings.PE1_Resolution*(settings.Matrix_Size(1)/2 - (settings.Matrix_Size(1)*(1-settings.PE1_Partial_Fourier)))); % (zero-indexed)
A = Centre_Line+1:1:(settings.Scan_Size(1)-1);
B = (Centre_Line)-1:-1:0;
N = min(numel(A),numel(B));
if strcmp(settings.PE1_Reordering,'CentricOut')
    ReorderPE1(1,:) = [Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE1(2,:) = ReorderPE1(1,:);
elseif strcmp(settings.PE1_Reordering,'CentricIn')
    ReorderPE1(1,:) = fliplr([Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1); % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE1(2,:) = ReorderPE1(1,:);
elseif strcmp(settings.PE1_Reordering,'CentricInOut')
    ReorderPE1(1,:) = fliplr([Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1); % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE1(2,:) = fliplr(ReorderPE1(1,:));
elseif strcmp(settings.PE1_Reordering,'LinearUp')
    ReorderPE1(1,:) = (0:1:round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1))-1) +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab)
    ReorderPE1(2,:) = ReorderPE1(1,:);
elseif strcmp(settings.PE1_Reordering,'LinearDown')
    ReorderPE1(1,:) = (0:1:round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1))-1) +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab)
    ReorderPE1(2,:) = ReorderPE1(1,:);
else
    error('Reordering scheme not recognised.');
end

% Now reorder PE2
if size(settings.Matrix_Size,2) > 1
    % Determine which line is the centre (i.e. not the centre of Matrix_Size for partial Fourier)
Centre_Line = ceil(settings.PE2_Resolution*(settings.Matrix_Size(2)/2 - (settings.Matrix_Size(2)*(1-settings.PE2_Partial_Fourier)))); % (zero-indexed)
A = Centre_Line+1:1:(round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2))-1);
B = (Centre_Line)-1:-1:0;
N = min(numel(A),numel(B));
if strcmp(settings.PE2_Reordering,'CentricOut')
    ReorderPE2(1,:) = [Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE2(2,:) = ReorderPE2(1,:);
elseif strcmp(settings.PE2_Reordering,'CentricIn')
    ReorderPE2(1,:) = fliplr([Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1); % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE2(2,:) = ReorderPE2(1,:);
elseif strcmp(settings.PE2_Reordering,'CentricInOut')
    ReorderPE2(1,:) = [Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    ReorderPE2(2,:) = fliplr(ReorderPE2(1,:));
elseif strcmp(settings.PE2_Reordering,'LinearUp')
    ReorderPE2(1,:) = (0:1:round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2))-1) +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab)
    ReorderPE2(2,:) = ReorderPE2(1,:);
elseif strcmp(settings.PE2_Reordering,'LinearDown')
    ReorderPE2(1,:) = (0:1:round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2))-1) +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab)
    ReorderPE2(2,:) = ReorderPE2(1,:);
else
    error('Reordering scheme not recognised.');
end
else
    ReorderPE2 = 'NA';
end


end

