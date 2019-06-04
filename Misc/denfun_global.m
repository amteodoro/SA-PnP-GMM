function [x_hate, PSNR_estimatee] = denfun_global(x, y, sigma, pd, ps, K,ld, sub, externalmix, supportvar)

y1 = y;

yy1 = wextend(2,'sym',y1,[pd-1,pd-1]);

[fN, fM] = size(y);

cont = 1;
PSNR_prev = 0;

count = 0;
while cont
    if sub == 0 || count == 1
        cont = 0;
    end
    
    sigma2 = sigma^2;
    
    global_dc = 0;
    
    yy = wextend(2,'sym',y,[pd-1,pd-1]);
    y_patches = im2colstep(yy,[pd,pd],[ps,ps]);
    y_patches_dc = mean(y_patches);
    y_patches_ac=bsxfun(@minus,y_patches,y_patches_dc);
    
    y_patches_train = y_patches_ac;
    [dim,num] = size(y_patches_ac);
    
    
    if ~isempty(ld)
        disp('Using mixture from database.')
        load(ld)
    else
        [prob,Scomp,Ucomp] = ...
            EM_zeromean(y_patches_train,K,sigma);
    end
    if ~exist('supportvar', 'var')
%         fprintf('No fixed support.\n')
        supportvar = [];
    end
    
    [x_hat_patches,post_cov] = GMM_inference(y_patches_ac,zeros(dim,K),prob,Scomp,Ucomp,sigma, supportvar);

    weights = (1./(post_cov+ eps));
    
    [x_haty, aux_best] = dccomp(x_hat_patches, y_patches_dc, weights, pd, ps, size(y), size(yy));
    
    
    PSNR_best = -10*log10(var( (x(:)-x_haty(:))));
        
    x_hatb = x_haty;
        
    if externalmix == 1
        p1 = prob;
        S1 = Scomp;
        U1 = Ucomp;
        
        l1 = length(p1);
                
        %% Try other databases
        database = 'text_trainpd8K20ps1';
        if exist(strcat(database, '.mat'), 'file')
            load(database);
        else
            database = 'text_train';
            extension = 'png';
            K = 20;
            [saveastext] = custommix(database, extension, pd, K, 1);
            load(saveastext);
        end
        
        textp = prob;
        textU = Ucomp;
        textS = Scomp;
        
        l2 = length(prob);
        
        database = 'mri_trainpd8K20ps4';
        
        if exist(strcat(database, '.mat'), 'file')
            load(database);
        else
            database = 'mri_train';
            extension = 'tif';
            K = 20;
            [saveastext] = custommix(database, extension, pd, K, 4);
            load(saveastext);
        end
        
        brainp = prob;
        brainU = Ucomp;
        brainS = Scomp;
        
        l3 = length(prob);
        
        database = 'algore_trainpd8K20ps1';
        if exist(strcat(database, '.mat'), 'file')
            load(database);
        else
            database = 'algore_train';
            extension = 'png';
            K = 20;
            [saveastext] = custommix(database, extension, pd, K, 1);
            load(saveastext);
        end
        
        facep = prob;
        faceU = Ucomp;
        faceS = Scomp;
        
        l4 = length(prob);
        
        database = 'fingerprints_trainpd8K20ps4';
        if exist(strcat(database, '.mat'), 'file')
            load(database);
        else
            database = 'fingerprints_train';
            extension = 'tif';
            K = 20;
            [saveastext] = custommix(database, extension, pd, K, 4);
            load(saveastext);
        end
        
        fingerp = prob;
        fingerU = Ucomp;
        fingerS = Scomp;
        
        l5 = length(prob);
        
        p = [p1; textp; brainp; facep; fingerp]/5;
        S = [S1, textS, brainS, faceS, fingerS];
        U = cat(3, U1, cat(3, textU, cat(3, brainU, cat(3, faceU, fingerU))));

        alpha_expansion = 0;
        
        for th = [0.3]
            
            for mod = 1:length(p)
                cov(:,:,mod) = U(:,:,mod)*bsxfun(@times,S(:,mod),U(:,:,mod)');
            end
            
            [indic] = computePost(y_patches_ac, p, zeros(dim, length(p)), cov + repmat(sigma^2*eye(pd^2),[1, 1, length(p)]));

            A = bsxfun(@minus,indic,max(indic,[],1));
            B = exp(A);
            normindic = bsxfun(@rdivide,B,sum(B,1));
            
            aux(1,:) = sum(normindic(1:l1,:));
            aux(2,:) = sum(normindic(l1+1:l2+l1,:));
            aux(3,:) = sum(normindic(l2+l1+1:l3+l2+l1,:));
            aux(4,:) = sum(normindic(l3+l2+l1+1:l4+l3+l2+l1,:));
            aux(5,:) = sum(normindic(l4+l3+l2+l1+1:end,:));
            
            if ~alpha_expansion
                
                [val,index] = max(aux);
                
                labels = ones(size(val));
                
                labels(val>th) = index(val>th);
                
            else
                aux_res = reshape(aux(1,:), [fN + pd - 1, fM + pd - 1]);
                aux_res(:,:,2) = reshape(aux(2,:), [fN + pd - 1, fM + pd - 1]);
                aux_res(:,:,3) = reshape(aux(3,:), [fN + pd - 1, fM + pd - 1]);
                aux_res(:,:,4) = reshape(aux(4,:), [fN + pd - 1, fM + pd - 1]);
                aux_res(:,:,5) = reshape(aux(5,:), [fN + pd - 1, fM + pd - 1]);
                
                probMatrix = log(aux_res + eps);
                beta = 1.5;
                Sc = ones(5) - eye(5);
                %Sc = [0 2 2 2 2; 2 0 1 1 1; 2 1 0 1 1; 2 1 1 0 1; 2 1 1 1 0];
                gch = GraphCut('open', -probMatrix, beta.*Sc);
                [gch, labels] = GraphCut('expand', gch);
                gch = GraphCut('close', gch);
                labels = labels + 1;
                
                
            end
            
%             imidx2 = reshape(labels, [fN + pd - 1, fM + pd - 1]);
%             figure, imagesc(imidx2)
%             drawnow
%             
            x_hat_patches = zeros(size(weights));
            
            if ~isempty(y_patches_ac(:,labels == 1))
                [x_hat_patches(:,labels == 1),post_cov(:,labels == 1)] = GMM_inference(y_patches_ac(:,labels == 1), zeros(dim,length(p1)), p1, S1, U1, sigma);
            end
            if ~isempty(y_patches_ac(:,labels == 2))
                [x_hat_patches(:,labels == 2),post_cov(:,labels == 2)] = GMM_inference(y_patches_ac(:,labels == 2), zeros(dim,length(textp)), textp, textS, textU, sigma);
            end
            if ~isempty(y_patches_ac(:,labels == 3))
                [x_hat_patches(:,labels == 3),post_cov(:,labels == 3)] = GMM_inference(y_patches_ac(:,labels == 3), zeros(dim,length(brainp)), brainp, brainS, brainU, sigma);
            end
            if ~isempty(y_patches_ac(:,labels == 4))
                [x_hat_patches(:,labels == 4),post_cov(:,labels == 4)] = GMM_inference(y_patches_ac(:,labels == 4), zeros(dim,length(facep)), facep, faceS, faceU, sigma);
            end
            if ~isempty(y_patches_ac(:,labels == 5))
                [x_hat_patches(:,labels == 5),post_cov(:,labels == 5)] = GMM_inference(y_patches_ac(:,labels == 5), zeros(dim,length(fingerp)), fingerp, fingerS, fingerU, sigma);
            end
            
            weights = 1./(post_cov+eps);
            
            [x_hat, aux_best] = dccomp(x_hat_patches, y_patches_dc, weights, pd, ps, size(y), size(yy));
            
            psnr = -10*log10(var( (x(:)-x_hat(:))));

            x_hatb = x_hat;
            PSNR_best = psnr;
            
        end
        
    end
    
    patch_vars = var(y_patches);
    PSNR_bestm = 0;
    
    for var_th = 0.5:0.02:1.5
        %         for var_th = 1
        
        aux_patches = x_hat_patches;
        weightsflat = weights;
        
        flat_patches = find(patch_vars < var_th*sigma2);
        
        aux_patches(:,flat_patches) = kron(y_patches_dc(flat_patches),ones(dim,1));
        
        %             weightsflat(:,flat_patches) = dim / sigma2;
        
        non_flat_patches = find(patch_vars >= var_th*sigma2);
        
        aux_patches(:,non_flat_patches) = bsxfun(@plus, aux_patches(:,non_flat_patches), y_patches_dc(non_flat_patches));
        
        aux_patches = min(max(aux_patches,-global_dc),1-global_dc);
        
        x_hat_raw = col2imstep(aux_patches.*(weightsflat),size(yy),[pd,pd],[ps,ps]);
        normalize = col2imstep((weightsflat),size(yy),[pd,pd],[ps,ps]);
        x_hat1 = x_hat_raw ./ normalize;
        
        x_hat = x_hat1(pd:(fN+pd-1),pd:(fM+pd-1)) + global_dc;
        
        x_hat = min(max(x_hat,0),1);
        
        PSNR_estim = -10*log10(var( (x(:)-x_hat(:))));
        
        if PSNR_estim > PSNR_best
            x_hatb = x_hat;
            aux_best = aux_patches;
            PSNR_best = PSNR_estim;
        end
        
        
    end
    x_hate = x_hatb;
    
    PSNR_estimatee = PSNR_best;
    
    
    x_hate1 = x_hate - global_dc;
    xxhate = wextend(2,'sym',x_hate1,[pd-1,pd-1]);
    

    patch_disagreement = aux_best - im2colstep(xxhate,[pd,pd],[ps,ps]);
    y_hat = im2colstep(yy1, [pd,pd], [ps,ps]) - patch_disagreement;
    %dis = col2imstep(patch_disagreement, size(yy), [pd,pd], [ps,ps]);
    y_hatraw = col2imstep(y_hat, size(yy), [pd,pd], [ps,ps]);
    normalize = col2imstep(ones(size(y_hat)),size(yy),[pd,pd],[ps,ps]);
    y_hat1 = y_hatraw ./ normalize;
    y = y_hat1(pd:(fN+pd-1),pd:(fM+pd-1)) + global_dc;
    
    if   PSNR_estimatee - PSNR_prev > 1e-2
        x_prev = x_hate;
        PSNR_prev = PSNR_estimatee;
    else
        cont = 0;
    end
    
    count = count + 1;
end

x_hate = x_prev;
PSNR_estimatee = PSNR_prev;



