from cbicall.goodbye import GoodBye


def test_goodbye_returns_known_word(monkeypatch):
    # Make the output deterministic
    monkeypatch.setattr("cbicall.goodbye.random.choice", lambda seq: seq[0])

    gb = GoodBye()
    word = gb.say_goodbye()
    assert word in GoodBye.WORDS
    assert word == GoodBye.WORDS[0]
