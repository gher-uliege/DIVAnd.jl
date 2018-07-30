
function sub2lin(sz,tmp)
    l = 0

    for i = length(tmp):-1:2
        l = (tmp[i]-1 + l)*sz[i-1]
    end

    l += tmp[1]
    return l
end

function lin2sub!(sz,ind::Number,sub)
    # tmp is here 0-based
    tmp = ind - 1

    for j = 1:length(sub)
        # integer division
        tmp2 = tmp ÷ sz[j]

        # make sub 1-based
        sub[j] = tmp - tmp2 * sz[j] +1
        tmp = tmp2
    end
end

sparse_diag(d) = sparse(1:length(d),1:length(d),d)

function sparse_pack(mask)
    j =
        if VERSION >= v"0.7.0-beta.0"
            LinearIndices(mask)[findall(mask)]
        else
            find(mask)
        end

    m = length(j)
    i = collect(1:m)
    s = ones(m)
    n = length(mask)
    H = sparse(i,j,s,m,n)
end
