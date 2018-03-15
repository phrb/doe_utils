using FileIO

results = RemoteChannel(() -> Channel(32))

for p in workers()
    @async remote_do(measure, p, results)
end

t = take!(results)

for i = 1:length(workers()) - 1
    t = merge(t, take!(results))
end

t

FileIO.save("./test.csv", t)
