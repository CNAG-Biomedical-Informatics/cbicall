import random


class GoodBye:
    """
    Just returns a random farewell, like the Perl GoodBye package.
    """

    WORDS = [
        "Aavjo",
        "Abar Dekha-Hobe",
        "Adeus",
        "Adios",
        "Aloha",
        "Alvida",
        "Ambera",
        "Annyong hi Kashipshio",
        "Arrivederci",
        "Auf Wiedersehen",
        "Au Revoir",
        "Ba'adan Mibinamet",
        "Dasvidania",
        "Donadagohvi",
        "Do Pobatchenya",
        "Do Widzenia",
        "Eyvallah",
        "Farvel",
        "Ha Det",
        "Hamba Kahle",
        "Hooroo",
        "Hwyl",
        "Kan Ga Waanaa",
        "Khuda Hafiz",
        "Kwa Heri",
        "La Revedere",
        "Le Hitra Ot",
        "Ma'as Salaam",
        "Mikonan",
        "Na-Shledanou",
        "Ni Sa Moce",
        "Paalam",
        "Rhonanai",
        "Sawatdi",
        "Sayonara",
        "Selavu",
        "Shalom",
        "Totsiens",
        "Tot Ziens",
        "Ukudigada",
        "Vale",
        "Zai Geen",
        "Zai Jian",
        "Zay Gesunt",
    ]

    def say_goodbye(self) -> str:
        return random.choice(self.WORDS)
