using BioSequences
using BioAlignments
using Statistics

function aln1(qd, td, so)
    problem = GlobalAlignment()
    s = AffineGapScoreModel(match=-1, mismatch=-2, gap_open=-4, gap_extend=-2)
    a = pairalign(problem, qd, td, s, score_only=(so != 0))
    if so == 0
        c = alignment(a)
    end
    1, score(a)
end

function main(query, target, so)
    fq = open(query)
    ft = open(target)
    total, score = 0, 0
    for (q, t) in zip(eachline(fq), eachline(ft))
        t, s = aln1(q, t, so)
        total += t
        score += s
    end
    close(fq)
    close(ft)
    total, score
end

function run()
    for so in 0:1
        et = []
        for i in 1:3
            e = @elapsed t, s = main(ARGS[1], ARGS[2], so)
            println("true $t $s $e")
            push!(et, e[0])
        end
        m = mean(et)
        s = std(et)
        println("[sw-time] julia $so $m $s")
    end
end

e = @elapsed run()
t = e[0]
println("total: $t")

