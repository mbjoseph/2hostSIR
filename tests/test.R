app <- ShinyDriver$new("..")
app$snapshotInit("test")

app$snapshot()
app$setInputs(b1 = 1)
app$snapshot()
app$setInputs(b1 = 10)
app$snapshot()
