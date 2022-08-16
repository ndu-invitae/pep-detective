import sqlite3

from pandas import read_sql_query

from pep_detective.src.ancova import AncovaResult


def init_conn(db: str) -> sqlite3.Connection:
    """create a database connection"""
    conn = None
    try:
        conn = sqlite3.connect(db)
    except sqlite3.Error as e:
        print(e)

    return conn


def create_db(conn: sqlite3.Connection):
    """create a database if not exists"""

    sql_create = """CREATE TABLE IF NOT EXISTS pep_stats (sample_id TEXT PRIMARYKEY, ph_covar CHAR, p_enhancer FLOAT, p_suppressor FLOAT, CONSTRAINT id_unique UNIQUE (sample_id)) """
    with conn:
        conn.execute(sql_create)


def insert_db(conn: sqlite3.Connection, ancova_results: AncovaResult):
    """insert anvoca results to created database; if sample id exists then replace it"""

    sql_insert = (
        """ INSERT OR REPLACE INTO pep_stats (sample_id, ph_covar, p_enhancer, p_suppressor) VALUES (?,?,?,?)"""
    )
    with conn:
        conn.execute(
            sql_insert,
            (
                ancova_results.sample_id,
                ancova_results.ph_covar,
                ancova_results.p_enhancer,
                ancova_results.p_suppressor,
            ),
        )


def get_db(conn):
    sql_get = """SELECT * from pep_stats """
    with conn:
        df_query = read_sql_query("SELECT * from pep_stats", conn)
        df_query["ph_covar"] = df_query["ph_covar"].apply(
            lambda x: False if x == b"\x00" else True
        )  # convert bit back to bool
        return df_query
