from database_init import CSFRNASourceDatabase


def create_schema(conn):
    db = CSFRNASourceDatabase()
    db.conn = conn
    db.create_database_schema()
