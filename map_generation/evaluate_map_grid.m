function B = evaluate_map_grid(model,N,E)

    r = [N(:),E(:),zeros(numel(N),1)];
    B = model.predict_map(r);
    B = reshape(B,size(N));

end